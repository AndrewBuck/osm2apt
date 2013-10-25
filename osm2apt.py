#!/usr/bin/env python

# Global variables that are user settable.
shoulderWidth = 1.0

import argparse
import math
import copy

from imposm.parser import OSMParser
from shapely.geometry import *
from shapely.geometry.polygon import orient
from shapely.ops import cascaded_union

overpassQueryFile = open('overpass_query.txt', 'w')
overpassQueryFile.write('data=\n\n[timeout:600];\n\n(\n')

def lookahead(iterable):
    it = iter(iterable)
    last = it.next()
    for val in it:
        yield last, False
        last = val
    yield last, True

def metersToDeg(meters):
    return (meters / (6371000 * 2 * 3.1415927)) * 360

def computeHeading(coord1, coord2):
    lat1 = coord1[1];  lon1 = coord1[0];
    lat2 = coord2[1];  lon2 = coord2[0];
    lat1 *= math.pi / 180.0
    lon1 *= math.pi / 180.0
    lat2 *= math.pi / 180.0
    lon2 *= math.pi / 180.0
    dLon = lon2 - lon1;
    y = math.sin(dLon) * math.cos(lat2)
    x = math.cos(lat1)*math.sin(lat2) - \
        math.sin(lat1)*math.cos(lat2)*math.cos(dLon)
    hdg = math.atan2(y, x) * 180.0 / math.pi

    if (hdg < 0):
        hdg+= 360.0

    return hdg

def headingToRunwayInt(heading):
    num = round(heading/ 10.0)
    if num == 0:
        num = 36

    return int(num)

def headingToRunwayString(heading):
    return '{0:02d}'.format(headingToRunwayInt(heading))

def normalizeHeading(heading):
    while heading < 0:
        heading += 360

    while heading > 360:
        heading -= 360

    return heading

def headingToDeg(heading):
    return (450 - heading) % 360.0

def computeTurnTo(pos1, origHeading, pos2):
    newHeading = computeHeading(pos1, pos2)
    return computeHeadingDelta(origHeading, newHeading)

def computeHeadingDelta(origHeading, newHeading):
    normalizeHeading(origHeading)
    normalizeHeading(newHeading)
    diff = math.fabs(origHeading - newHeading)
    if diff <= 180:
        amount = diff
        if newHeading >= origHeading:
            direction = 'right'
        else:
            direction = 'left'
    else:
        amount = 360 - diff
        if origHeading >= newHeading:
            direction = 'right'
        else:
            direction = 'left'

    return (amount, direction)

def computeSegmentHeading(node, nodes, coords):

    # Determine which node along the way the node occurs on.
    index = nodes.index(node)

    # If the node is any but the last one, we use the direction of the
    # _following_ segment of the taxiway to get the direction, however
    # if the node is the last in the taxiway, then there is no
    # following segement so we use the previous segment instead.
    if index < (len(nodes) - 1):
        p1 = coords[index]
        p2 = coords[index+1]
        pos = p1
    else:
        p1 = coords[index-1]
        p2 = coords[index]
        pos = p2

    heading = computeHeading(p1, p2)

    return (heading, pos)

def computeHeadingAtPoint(line, point, towardPoint=-1, delta=metersToDeg(0.5)):
    distance = line.project(point)

    if towardPoint != -1:
        distanceToward = line.project(towardPoint)
        if distanceToward < distance:
            delta *= -1.0
            
    deltaPoint = line.interpolate(distance+delta)

    return computeHeading((point.x, point.y), (deltaPoint.x, deltaPoint.y))

def computeNearestObject(obj, otherObjs):
    tempDistance = -1
    shortestDistance = -1
    nearestObject = -1
    if isinstance(obj, SpatialObject):
        geometry = obj.geometry
    else:
        geometry = obj

    for otherObj in otherObjs:
        if isinstance(otherObj, SpatialObject):
            otherGeometry = otherObj.geometry
        else:
            otherGeometry = otherObj

        tempDistance = geometry.distance(otherGeometry)
        if shortestDistance == -1:
            shortestDistance = tempDistance
            nearestObject = otherObj
        else:
            if tempDistance < shortestDistance:
                shortestDistance = tempDistance
                nearestObject = otherObj

    return (nearestObject, shortestDistance)

def computeJunctionSigns(coord, ways):
    signs = []
    junctions = []
    wayWidths = []
    wayGeoms = []
    # Make lists of all the taxiway widths and centerline geometries so we can
    # find the maximum width and combine all the geometries into a temporary
    # union of all the taxiways/runways that join at this junction.
    for way in ways:
        wayWidths.append(way.width)
        wayGeoms.append(way.geometry)

    geometryUnion = cascaded_union(wayGeoms)

    # Compute the intersection points of the geometries, first at a large ring
    # out to where the signs will be placed, then again at a 1 meter ring right
    # around the junction node to see what direction each way enters the
    # junction from.
    setbackDistance = metersToDeg(max(wayWidths) + 5.0)
    setbackRing = Point(coord).buffer(setbackDistance).exterior
    setbackPoints = setbackRing.intersection(geometryUnion)

    junctionDistance = metersToDeg(1.0)
    junctionRing = Point(coord).buffer(junctionDistance).exterior
    junctionPoints = junctionRing.intersection(geometryUnion)
    
    # If at the setback distance there is only a single point rather than a
    # multipoint then at most one taxiway is long enough to even reach the
    # setback ring, so we just return no signs since they are not needed
    # anyway.
    if not isinstance(setbackPoints, MultiPoint):
        return signs

    for setbackPoint in setbackPoints:
        closestWay, distance = computeNearestObject(setbackPoint, ways)
        wayDistanceSetback = closestWay.geometry.project(setbackPoint)
        wayDistanceJunction = closestWay.geometry.project(Point(coord))
        if wayDistanceSetback > wayDistanceJunction:
            junctionPoint = closestWay.geometry.interpolate(wayDistanceJunction + metersToDeg(1.0))
        else:
            junctionPoint = closestWay.geometry.interpolate(wayDistanceJunction - metersToDeg(1.0))
        #junctionPoint, distance = computeNearestObject(closestWay, junctionPoints)
        setbackHeading = computeHeadingAtPoint(closestWay.geometry, setbackPoint, junctionPoint)
        junctionHeading = computeHeadingAtPoint(closestWay.geometry, junctionPoint, Point(coord))
        junctions.append((closestWay, setbackPoint, setbackHeading, junctionHeading))

    for (way1, setbackPoint1, setbackHeading1, junctionHeading1) in junctions:
        signLoc = travel(setbackPoint1.coords[0], setbackHeading1-90, metersToDeg(way1.width/2.0 + 2.5))

        subsignParts = []
        # Determine the text to place on the sign.
        for (way2, setbackPoint2, setbackHeading2, junctionHeading2) in junctions:
            if setbackPoint1.distance(setbackPoint2) > metersToDeg(1.0):
                headingString = ''
                deltaHeading, direction = computeHeadingDelta(junctionHeading1, junctionHeading2+180)
                directionLetter = direction[0]

                if deltaHeading <= 22.5:
                    headingString = '{^u}'
                elif deltaHeading <= 67.5:
                    headingString = '{^' + directionLetter + 'u}'
                elif deltaHeading <= 112.5:
                    headingString = '{^' + directionLetter + '}'
                elif deltaHeading <= 157.5:
                    headingString = '{^' + directionLetter + 'd}'
                else:
                    headingString = '{^d}'

                # Print leftward arrows to the left of the name and rightward arrows on the right.
                if direction == 'left':
                    deltaHeading *= -1.0
                    text = headingString + way2.name
                else:
                    text = way2.name + headingString

                subsignParts.append((deltaHeading, text, way2))

        text = ''
        # TODO: First sort subsignParts from lowest to highest deltaHeading.
        for (deltaHeading, partText, way), isLast in lookahead(subsignParts):
            typeString = ''
            if isinstance(way, Taxiway):
                typeString = '{@Y}'
            elif isinstance(way, Runway):
                typeString = '{@R}'

            text += typeString + partText

            if not isLast:
                text += '|'
        signs.append(Sign(signLoc, normalizeHeading(setbackHeading1), 2, text))

    return signs

def travel(startPos, heading, distance):
    heading = normalizeHeading(heading)
    x = startPos[0]
    y = startPos[1]
    angle = headingToDeg(heading) * math.pi / 180.0
    x += distance * math.cos(angle)
    y += distance * math.sin(angle)

    return (x, y)

def addOverpassQuery(type, id):
    if type == 'node':
        overpassQueryFile.write('  node({0});\n'.format(id))
    elif type == 'way':
        overpassQueryFile.write('  way({0});    >;\n'.format(id))
    elif type == 'relation':
        overpassQueryFile.write('  relation({0});    >;\n'.format(id))

def checkNodesFormClosedWay(nodes):
    if nodes[0] != nodes[-1]:
        print 'WARNING: Way not closed.  First node: %s  Last node: %s' % (nodes[0], nodes[-1])
        return False
    else:
        return True

def nodesToCoords(nodes, coordDict):
    coords = []
    for n in nodes:
        coords.append((coordDict[n][0], coordDict[n][1], n))

    return coords

def printLine(line, lineType, lineName):
    ret = ''

    if isinstance(line, MultiLineString) or isinstance(line, GeometryCollection):
        for c in line.geoms:
            ret += printLine(c, lineType, lineName)

        return ret

    if isinstance(line, LineString):
        ret += printLineSegment(line, lineType, lineName)
    else:
        print 'ERROR: Should never get here. (printLine())'

    return ret

def printLineSegment(line, lineType, lineName):
    ret = '120 {0}\n'.format(lineName)

    for coord, isLast in lookahead(line.coords):
        if not isLast:
            ret += '111 {0} {1} {2}\n'.format(coord[1], coord[0], lineType)
        else:
            ret += '115 {0} {1}\n'.format(coord[1], coord[0])

    return ret

def printArea(area, areaType, surface, smoothness, heading, name):
    ret = ''

    # TODO: Fix this to actually handle these if necessary, for now any time
    # this gets passed a GeometryCollection the thing is empty anyway so
    # nothing is lost by just exiting as we do now.
    if isinstance(area, GeometryCollection) or isinstance(area, MultiPolygon):
        for obj in area.geoms:
            ret += printArea(obj, areaType, surface, smoothness, heading, name)
        return ret

    if isinstance(area, Polygon):
        ret += printPolygon(area, areaType, surface, smoothness, heading, name)

    return ret

def printPolygon(area, areaType, surface, smoothness, heading, name):
    ret = '{0} {1} {2} {3} {4}\n'.format(areaType, surface, smoothness, heading, name)
    if area.exterior.is_ccw:
        coords = area.exterior.coords
    else:
        coords = reversed(area.exterior.coords)

    for coord, isLast in lookahead(coords):
        if not isLast:
            ret += '111 {0} {1}\n'.format(coord[1], coord[0])
        else:
            ret += '113 {0} {1}\n'.format(coord[1], coord[0])

    for ring in area.interiors:
        if ring.is_ccw:
            coords = reversed(ring.coords)
        else:
            coords = ring.coords

        for coord, isLast in lookahead(coords):
            if not isLast:
                ret += '111 {0} {1}\n'.format(coord[1], coord[0])
            else:
                ret += '113 {0} {1}\n'.format(coord[1], coord[0])

    return ret

def cut(line, distance):
    # Cuts a line in two at a distance from its starting point
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [
                LineString(coords[:i+1]),
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]),
                LineString([(cp.x, cp.y)] + coords[i:])]

def isAbandoned(tags):
    if 'abandoned' in tags:
        if tags['abandoned'] == 'yes':
            return True

    return False

def convertToUnit(string, unit):
    conversionFactors = {
        'ft' : { 'ft' : 1.0,
                 'm'  : 3.28084,
                 ''   : 3.28084},

        'm'  : { 'm'  : 1.0,
                 'ft' : 0.3048,
                 ''   : 1.0}
    }

    originalString = string
    string = string.strip()

    stringParts = string.split(None, 1)
    length = len(stringParts)

    if length == 0:
        print 'ERROR: Parsing units on string "%s" failed, returning 0.' % originalString
        return 0
    elif length == 1:
        stringParts.append('')

    stringNumber = float(stringParts[0])
    stringUnit = stringParts[1]

    if unit in conversionFactors:
        conversionRatioDict = conversionFactors[unit]
    else:
        print 'ERROR: Resultant unit of "%s" not found in conversion factors, returning 0.' % unit
        return 0

    if stringUnit in conversionRatioDict:
        conversionRatio = float(conversionRatioDict[stringUnit])
    else:
        print 'ERROR: Source unit of "%s" not found in conversion factors, returning 0.' % stringUnit
        print 'Acceptable source units and corresponding conversion ratios to target unit "%s" are:' % unit
        print conversionRatioDict
        return 0

    return float(stringNumber * conversionRatio)

def surfaceStringToInt(st):

    surfaceTypeDict = {'asphalt' : 1,
                       'paved' : 1,
                       'concrete' : 2,
                       'turf' : 3,
                       'grass' : 3,
                       'unpaved' : 3,
                       'dirt' : 4,
                       'gravel' : 5,
                       'dry_lakebed' : 12,
                       'water' : 13,
                       'snow' : 14,
                       'transparent' : 15}

    if st in surfaceTypeDict:
        return surfaceTypeDict[st]
    else:
        return surfaceTypeDict['unpaved']

# Takes a tags dictionary and a list of keys and returns the value of the first
# key in the dict that matches one from the list.  If no key is found in the
# list of keys given, then the value passed to 'default' is returned instead.
def coalesceValue(keys, tags, default):
    for key in keys:
        if key in tags:
            return tags[key]

    return default

class SpatialObject(object):
    
    def __init__(self):
        self.geometry = null

    def buildGeometry(self, coordDict, nodeDict):
        return

class AerodromeObject(SpatialObject):

    def __init__(self):
        self.parentAerodrome = null

class Aerodrome(SpatialObject):

    def __init__(self, name, code, ele, nodes):
        self.name = name
        self.code = code
        self.ele = convertToUnit(ele, 'ft')
        self.nodes = nodes
        self.assosciatedObjects = []
        # TODO: set viewpoint location, row code 14, one per aerodrome.

    def buildGeometry(self, coordDict, nodeDict):
        if checkNodesFormClosedWay(self.nodes):
            self.geometry = Polygon(nodesToCoords(self.nodes, coordDict))
        else:
            print 'Not building geometry for aerodrome since it is not closed.'

    def listObjectsByType(self, objType):
        ls = []
        for obj in self.assosciatedObjects:
            if isinstance(obj, objType):
                ls.append(obj)

        return ls

    def cleanGeometries(self):
        self.taxiwaySurfaces = {}
        for taxiway in self.listObjectsByType(Taxiway):
            for apron in self.listObjectsByType(Apron):
                # Cut away the parts of the taxiways surface as well as just the taxiway edge lines that are on aprons.
                taxiway.concreteGeometry = taxiway.concreteGeometry.difference(apron.geometry)
                taxiway.leftDotted = apron.geometry.intersection(taxiway.leftDotted)
                taxiway.rightDotted = apron.geometry.intersection(taxiway.rightDotted)
                taxiway.leftEdgeLine = taxiway.leftEdgeLine.difference(apron.geometry)
                taxiway.rightEdgeLine = taxiway.rightEdgeLine.difference(apron.geometry)


            for runway in self.listObjectsByType(Runway):
                # Cut away the parts of the taxiways surface as well as just the taxiway edge lines that are on runways.
                taxiway.leftEdgeLine = taxiway.leftEdgeLine.difference(runway.geometryPolygon)
                taxiway.rightEdgeLine = taxiway.rightEdgeLine.difference(runway.geometryPolygon)
                #taxiway.concreteGeometry = taxiway.concreteGeometry.difference(runway.geometryPolygon)

            # Combine all of the taxiways of the same surface type together into a single multipolygon.
            # TODO: This could probably be re-factored to use a 'cascading union' from shapely
            if taxiway.surfaceInteger in self.taxiwaySurfaces:
                self.taxiwaySurfaces[taxiway.surfaceInteger] = self.taxiwaySurfaces[taxiway.surfaceInteger].union(taxiway.concreteGeometry)
            else:
                self.taxiwaySurfaces[taxiway.surfaceInteger] = taxiway.concreteGeometry

            # Loop over each pair of taxiways and use a buffer along the center
            # of each taxiway to trim out the shoulder markings on the other
            # one; this cleans up the intersections.
            for taxiway2 in self.listObjectsByType(Taxiway):
                if taxiway2 is not taxiway:
                    center = taxiway.geometry.buffer(metersToDeg(taxiway.width / 2.0 - shoulderWidth))
                    taxiway2.leftEdgeLine = taxiway2.leftEdgeLine.difference(center)
                    taxiway2.rightEdgeLine = taxiway2.rightEdgeLine.difference(center)

    def toString(self):
        # Print out the main airport header line
        tempString = '1 {0} 0 0 {1} {2}\n'.format(self.ele, self.code, self.name)

        # Print out the boundary area of the aerodrome.
        tempString += printArea(self.geometry, 130, '', '', '', self.name)

        # Print out all of the assosciated objects for this airport to apt.dat.
        taxiwayCoords = {}
        for obj in self.assosciatedObjects:
            tempString += obj.toString()

        # Print out the combined taxiway surfaces.
        for surf, area in self.taxiwaySurfaces.items():
            area = area.buffer(0)
            area = area.simplify(metersToDeg(0.5))
            tempString += printArea(area, 110, surf, 0.15, 360, 'taxiway surface ' + str(surf))

        # If the airport has associated taxiways, also print out a taxiway network.
        taxiways = self.listObjectsByType(Taxiway)
        runways = self.listObjectsByType(Runway)
        taxiwaysAndRunways = taxiways[:] + runways[:]
        if len(taxiwaysAndRunways) > 0:
            # Start by building a set of the nodes used for taxiways at this specific airport.
            for taxiway in taxiwaysAndRunways:
                for coord in taxiway.coords:
                    if coord in taxiwayCoords:
                        taxiwayCoords[coord].append(taxiway)
                    else:
                        taxiwayCoords[coord] = [taxiway]

            # Now print out the lines for all of the taxiway/runway nodes for this airport.
            tempString += '1200\n'
            for coord in taxiwayCoords:
                tempString += '1201 {0} {1} both {2}\n'.format(coord[1], coord[0], coord[2])

            # Finally, print out the actual edges between the nodes.
            for way in taxiwaysAndRunways:
                prevCoord = 0
                for coord in way.coords:
                    if prevCoord == 0:
                        prevCoord = coord
                    else:
                        # TODO: Implement oneway taxiways.
                        if isinstance(way, Runway):
                            typeString = 'runway'
                        else:
                            typeString = 'taxiway'

                        tempString += '1202 {0} {1} twoway {2} {3}\n'.format(prevCoord[2], coord[2], typeString, way.name)

                        if typeString == 'runway':
                            tempString += '1204 arrival {0}\n'.format(way.runwayEndNames[0] + ',' + way.runwayEndNames[1])

                        prevCoord = coord

            # Place signs at the taxiway intersections
            for coord, ways in taxiwayCoords.items():
                if len(ways) >= 2:
                    signs = computeJunctionSigns(coord, ways)
                    for sign in signs:
                        tempString += sign.toString()

        return tempString

class Runway(AerodromeObject):

    def __init__(self, name, runwayEndNames, runwayEndNodes, surface, width, nodes):
        self.name = name
        self.runwayEndNames = runwayEndNames
        self.runwayEndNodes = runwayEndNodes
        self.surface = surface
        self.surfaceInteger = surfaceStringToInt(self.surface)
        self.width = convertToUnit(width, 'm')
        self.nodes = nodes

    def buildGeometry(self, coordDict, nodeDict):
        self.geometry = LineString(nodesToCoords(self.nodes, coordDict))
        self.coords = nodesToCoords(self.nodes, coordDict)

        # Check to see if the first end of the runway is drawn with the lower
        # numbered end first or if it needs to be reversed so that the
        # direction of the line lines up with the name.  Note that this assumes
        # the ref tag gives the lower number first, like 18/36.
        firstEnd = headingToRunwayInt(computeHeading(self.geometry.coords[0], self.geometry.coords[-1]))
        if firstEnd > 18:
            self.geometry.coords = list(self.geometry.coords)[::-1]
            self.runwayEndNodes = list(self.runwayEndNodes)[::-1]

        self.heading = computeHeading(self.geometry.coords[0], self.geometry.coords[-1])
        self.geometryPolygon = self.geometry.buffer(metersToDeg(self.width/2.0), 1)

    def toString(self):
        return '100 {0} {1} 0 0.15 0 0 1 {2} {3} {4} 0 0 1 0 0 0 {5} {6} {7} 0 0 1 0 0 0\n'.format(self.width, self.surfaceInteger, self.runwayEndNames[0], self.geometry.coords[0][1], self.geometry.coords[0][0], self.runwayEndNames[1], self.geometry.coords[-1][1], self.geometry.coords[-1][0])

class Taxiway(AerodromeObject):

    def __init__(self, name, surface, width, nodes):
        self.name = name
        self.surface=surface
        self.surfaceInteger = surfaceStringToInt(self.surface)
        self.width = convertToUnit(width, 'm')
        self.nodes = nodes
        self.holdingPositions = []

    def buildGeometry(self, coordDict, nodeDict):
        self.coords = nodesToCoords(self.nodes, coordDict)
        self.geometry = LineString(self.coords)
        self.leftEdgeLine = self.geometry.parallel_offset(metersToDeg(self.width / 2.0 - shoulderWidth), 'left', 2)
        self.rightEdgeLine = self.geometry.parallel_offset(metersToDeg(self.width / 2.0 - shoulderWidth), 'right', 2)
        self.leftDotted = copy.deepcopy(self.leftEdgeLine)
        self.rightDotted = copy.deepcopy(self.rightEdgeLine)
        self.concreteGeometry = orient(self.geometry.buffer(metersToDeg(self.width)/2.0, 2))

        # Loop over all the nodes that make up the taxiway and make a list of
        # any that are holding positions.
        # TODO: Currently only checks normal holds, should check for ILS holds as well.
        for node in self.nodes:
            if node in nodeDict:
                nodeProps = nodeDict[node]
                if 'aeroway' in nodeProps[0]:
                    if nodeProps[0]['aeroway'] == 'holding_position':
                        self.holdingPositions.append((node, nodeProps[0], nodeProps[1]))

    def toString(self):
        ret = ''

        # Print out holding position lines on top of the concrete (these lines
        # are purely visual, X-plane does not respect these in taxi routing)
        self.signs = []
        for node, tags, coord in self.holdingPositions:
            headingAndPosition = computeSegmentHeading(node, self.nodes, self.coords)
            hdg = headingAndPosition[0]
            pos = headingAndPosition[1]
            # TODO: Need to properly set runway.
            ret += Sign(travel(pos, hdg+90, metersToDeg(self.width/2.0 + 2.5)), hdg, 3, '{@R}18-36').toString()
            ret += Sign(travel(pos, hdg-90, metersToDeg(self.width/2.0 + 2.5)), hdg, 3, '{@R}18-36').toString()

            # TODO: The dashed yellow lines should be on the left side of this line, which should be closer to the runway.  Currently they are drawn arbitrarily on one side, might need to reverse start and end if they are on the wrong side for a particular runway.
            lineStart = travel(pos, hdg-90, metersToDeg(self.width/2.0 - shoulderWidth))
            lineEnd = travel(pos, hdg+90, metersToDeg(self.width/2.0 - shoulderWidth))

            ret += printLine(LineString((lineStart, lineEnd)), 4, 'hold line {0}\n'.format(self.name))

        # Print out taxiway centerlines and shoulder lines on top of the concrete.
        ret += printLine(self.geometry, 1, 'taxiway {0} centerline'.format(self.name))
        ret += printLine(self.leftEdgeLine, 3, 'taxiway {0} left line'.format(self.name))
        ret += printLine(self.rightEdgeLine, 3, 'taxiway {0} right line'.format(self.name))
        ret += printLine(self.leftDotted, 2, 'taxiway {0} dotted line'.format(self.name))
        ret += printLine(self.rightDotted, 2, 'taxiway {0} dotted line'.format(self.name))
        return ret

class Windsock(AerodromeObject):

    def __init__(self, coord, lit):
        self.coord = coord
        if lit in ('yes', 'true', '1'):
            self.lit = 1
        else:
            self.lit = 0

    def buildGeometry(self, coordDict, nodeDict):
        self.geometry = Point(self.coord)

    def toString(self):
        return '19 {0} {1} {2} WS\n'.format(self.coord[1], self.coord[0], self.lit)

class Apron(AerodromeObject):

    def __init__(self, name, nodes, surface):
        self.name = name
        self.nodes = nodes
        self.surface = surface
        self.surfaceInteger = surfaceStringToInt(self.surface)

    def buildGeometry(self, coordDict, nodeDict):
        if checkNodesFormClosedWay(self.nodes):
            # TODO: Need to make sure the coords in the resulting geometry for a counter clockwise ring.
            self.geometry = orient(Polygon(nodesToCoords(self.nodes, coordDict)))
        else:
            print 'Not building geometry for apron since it is not closed.'

    def toString(self):
        ret = printArea(self.geometry, 110, self.surfaceInteger, 0.15, 360, self.name)

        # If this apron is named, also create a start location for it.
        # TODO: If an airport has no named aprons we should choose the largest one and put a start location there.  An alternative is for the aerodrome class to make sure that it has at least one named apron by naming the largest on 'Ramp' or something if none are named.
        if len(self.name) > 0:
            startupLoc = self.geometry.representative_point()
            # TODO: Set the startup heading to be something reasonable, toward the nearest taxiway would be a good one.
            ret += '15 {0} {1} 360 {2}\n'.format(startupLoc.y, startupLoc.x, self.name)

        return ret

class Beacon(AerodromeObject):

    def __init__(self, coord, color):
        self.coord = coord
        self.color = color

    def buildGeometry(self, coordDict, nodeDict):
        self.geometry = Point(self.coord)

    def toString(self):
        return '18 {0} {1} {2} BCN\n'.format(self.coord[1], self.coord[0], self.color)

class ParkingPosition(SpatialObject):

    def __init__(self, coord, positionTypeTag, width):
        self.coord = coord
        self.positionTypeTag = positionTypeTag
        self.heading = 360
        self.typeString = 'tie-down'
        self.aircraftTypeString = 'props'
        self.width = width

        if self.positionTypeTag == 'parking_position':
            self.typeString = 'tie_down'
            self.aircraftTypeString = 'props'
        elif self.positionTypeTag == 'hangar':
            self.typeString = 'hangar'
            self.aircraftTypeString = 'all'
        elif self.positionTypeTag == 'gate':
            self.typeString = 'gate'
            self.aircraftTypeString = 'jets|heavy'
        else:
            print 'ERROR: Parking position type "%s" unknown, may be exported incorrectly.'

    def buildGeometry(self, coordDict, nodeDict):
        self.geometry = Point(self.coord)
        self.leftWingTip = travel(self.coord, self.heading-90, metersToDeg(self.width/2.0))
        self.rightWingTip = travel(self.coord, self.heading+90, metersToDeg(self.width/2.0))
        self.tail = travel(self.coord, self.heading+180, metersToDeg(self.width/2.0))
        self.wingLine = LineString((self.leftWingTip, self.rightWingTip))
        self.tailLine = LineString((self.coord, self.tail))

    def toString(self):
        ret = '1300 {0} {1} {2} {3} {4}\n'.format(self.coord[1], self.coord[0], round(self.heading), self.typeString, self.aircraftTypeString)
        ret += printLine(self.wingLine, 1, self.typeString)
        ret += printLine(self.tailLine, 1, self.typeString)
        return ret

class LightedObject(AerodromeObject):

    def __init__(self, coord, typeName):
        self.coord = coord
        self.typeName = typeName
        self.heading = 360
        # TODO: Need to figure out tag for glideslope and read it if present.
        self.glideslope = 3.0
        self.runway = '01'

    def buildGeometry(self, coordDict, nodeDict):
        self.geometry = Point(self.coord)

    def toString(self):
        # Determine the nearest runway to this object.
        nearestRunway, shortestDistance = computeNearestObject(self.geometry, self.parentAerodrome.listObjectsByType(Runway))
        turn = -1
        if nearestRunway != -1:
            # Find out what distance along the runway this object is at.
            distance = nearestRunway.geometry.project(self.geometry)
            runwayFraction = distance / nearestRunway.geometry.length
            runwayLoc = nearestRunway.geometry.interpolate(distance)
            # Determine the true heading of the runway and also work out which
            # end of the runway the light is on.
            if runwayFraction < 0.5:
                self.heading = computeHeading(nearestRunway.geometry.coords[0], nearestRunway.geometry.coords[-1])
                self.runway = nearestRunway.runwayEndNames[0]
            else:
                self.heading = computeHeading(nearestRunway.geometry.coords[-1], nearestRunway.geometry.coords[0])
                self.runway = nearestRunway.runwayEndNames[1]

            self.heading = int(round(self.heading))
            if self.heading == 0:
                self.heading = 360

            # If we are moving down the runway and are at the location with the
            # light directly beside us, compute which way (left or right) we
            # would need to turn to face the light, i.e. determine which side
            # of the runway it is on.
            turn = computeTurnTo(runwayLoc.coords[:][0], self.heading, self.geometry.coords[:][0])

        # Determine the appropriate code to use for X-plane to indicate the type of the light.
        if self.typeName == 'vasi':
            self.typeCode = 1
        elif self.typeName == 'papi':
            if turn != -1 and turn[1] == 'left':
                self.typeCode = 2
            else:
                self.typeCode = 3

        else:
            print 'ERROR: Lighted object of type "%s" is unknown, skipping output of this object.' % self.typeName

        return '21 {0} {1} {2} {3} {4} {5} {6}\n'.format(self.coord[1], self.coord[0], self.typeCode, round(self.heading), self.glideslope, self.runway, self.typeName)

class Sign(AerodromeObject):

    def __init__(self, coord, heading, size, text):
        self.coord = coord
        self.heading = heading
        self.size = size
        self.text = text

    def buildGeometry(self):
        self.geometry = Point(self.coord)

    def toString(self):
        return '20 {0} {1} {2} 0 {3} {4}\n'.format(self.coord[1], self.coord[0], round(self.heading), self.size, self.text)

class Osm2apt_class(object):

    aerodromes = []

    # Every one of the lists in this next section should be added to
    # objectLists like 'runways' is.
    objectLists = []
    runways = [];    objectLists.append(runways)
    taxiways = [];   objectLists.append(taxiways)
    windsocks = [];   objectLists.append(windsocks)
    aprons = [];   objectLists.append(aprons)
    beacons = [];   objectLists.append(beacons)
    lightedObjects = [];   objectLists.append(lightedObjects)
    parkingPositions = [];   objectLists.append(parkingPositions)

    coords = []
    coordDict = {}
    nodes = []
    nodeDict = {}

    # Callback method to simply read in all node ID's and coordinates
    def coordsCallback(self, coords):
        for c in coords:
            self.coords.append(c)

    # Callback method to read all nodes which have tags on them
    def nodesCallback(self, nodes):
        for osmid, tags, coord in nodes:
            self.nodes.append((osmid, tags, coord))

            if 'aeroway' in tags:
                #node: aeroway=windsock
                if tags['aeroway'] == 'windsock':
                    lit = coalesceValue(('lit'), tags, 'no')
                    self.windsocks.append(Windsock(coord, lit))

                #node: aeroway=vasi|papi
                if tags['aeroway'] == 'vasi' or tags['aeroway'] == 'papi':
                    self.lightedObjects.append(LightedObject(coord, tags['aeroway']))

                #node: aeroway=parking_position
                if tags['aeroway'] == 'parking_position' \
                or tags['aeroway'] == 'hangar' \
                or tags['aeroway'] == 'gate':
                    width = convertToUnit(coalesceValue(('width'), tags, '10 m'), 'm')
                    self.parkingPositions.append(ParkingPosition(coord, tags['aeroway'], width))

            if 'man_made' in tags:
                #node: man_made=beacon
                if tags['man_made'] == 'beacon':
                    # TODO: Need to work out the tagging (if it even exists) for the beacon color type.
                    color = 1
                    self.beacons.append(Beacon(coord, color))

    # Callback method to process ways
    def waysCallback(self, ways):

        for osmid, tags, refs in ways:

            if 'aeroway' in tags:
                # way: aeroway=aerodrome
                if tags['aeroway'] == 'aerodrome':
                    print '\nFound an aerodrome.  Tags:\n', tags

                    # Get the airport name, if no name is listed use the
                    # default name listed below, no need to skip it in this case.
                    if 'name' in tags:
                        aerodromeName = tags['name']
                    else:
                        aerodromeName = 'Unnamed Airport'

                    print 'Airport name is: ', aerodromeName

                    # Get the ICAO code first, if that does not exist fall back
                    # on the FAA code, if neither exist we cannot import the
                    # airport so just skip it.
                    aerodromeCode = coalesceValue(('icao', 'faa', 'ref'), tags, '')
                    if aerodromeCode == '':
                        print 'ERROR: Aerodrome does not have an ICAO or FAA code, skipping.  Way ID: ', osmid
                        addOverpassQuery('way', osmid)
                        continue

                    print 'Airport Code is: ', aerodromeCode

                    # Get the airport elevation, if there is no 'ele' tag skip
                    # the airport.
                    if 'ele' in tags:
                        aerodromeELE = tags['ele']
                    else:
                        print 'ERROR: Aerodrome does not have an elevation (ele), skipping. Way ID: ', osmid
                        addOverpassQuery('way', osmid)
                        continue

                    print 'Airport Elevation is: ', aerodromeELE

                    # We have successfully read all the data for this aerodrome
                    # so add it to the list of completed aerodromes to be put
                    # in the output file.
                    self.aerodromes.append(Aerodrome(aerodromeName, aerodromeCode, aerodromeELE, refs))

                # way: aeroway=runway
                elif tags['aeroway'] == 'runway':
                    print '\nFound a runway.  Tags:\n', tags

                    # Check to make sure the runway is not a closed loop, that
                    # is, make sure the first and last nodes are not the same node.
                    if (refs[0] != refs[-1]):
                        runwayEndNodes = [refs[0], refs[-1]]
                    else:
                        print 'ERROR: Runway is a closed loop, skipping.  Way ID: ', osmid
                        addOverpassQuery('way', osmid)
                        continue

                    # Check to see if the runway is abandonded, if so skip it.
                    if isAbandoned(tags):
                        print 'NOTICE: Runway is abandoned, skipping. Way ID: ', osmid
                        continue

                    # Use the 'ref' tag to get the runway end names (for most
                    # runways this should be filled in).  We also store the
                    # first and last node ID at this point, although these may
                    # need to be reversed later to make the actual coords match
                    # up with the end names.
                    # TODO: Right now we skip the runway if this is not set, in
                    # the future we could guess the runway name based on its orientation.
                    name = coalesceValue(('ref', 'name'), tags, 'unnamed')
                    if name != 'unnamed':
                        # TODO: Need to check to see what the runway separator is, most use '/' but some runways use a '-' and others are probably used as well.
                        runwayEndNames = name.split('/')
                        print 'Runway end names are: ', runwayEndNames
                    else:
                        print 'ERROR: Runway does not have a "ref" or "name" tag set, skipping.  Way ID: ', osmid
                        addOverpassQuery('way', osmid)
                        continue

                    # Determine the runway surface type or assume 'concrete' if not specified.
                    if 'surface' in tags:
                        runwaySurface = tags['surface']
                    else:
                        runwaySurface = ''

                    print 'Runway surface is: ', runwaySurface

                    # Determine the runway width or use a sensible default.
                    if 'width' in tags:
                        runwayWidth = tags['width']
                    else:
                        runwayWidth = '20 m'

                    print 'Runway width is: ', runwayWidth

                    #TODO: Process the 'length' field and set appropriately, currently
                    # the length will be based on the coords of the end nodes on the way.

                    #TODO: Determine and set the runway shoulder type, not sure
                    # if this is even tagged in OSM yet.

                    #TODO: Determine and set the runway smoothness (not used in X-plane yet).

                    #TODO: Determine and set the runway centerline lighting.

                    #TODO: Determine and set the runway edge lighting.

                    #TODO: Handle displaced thresholds and blast pads.

                    #TODO: Handle runway markings, TDZ lighting, and REIL lighting.

                    # We have successfully read all the data for this runway
                    # so add it to the list of completed runways to be put
                    # in the output file.
                    self.runways.append(Runway(name, runwayEndNames, runwayEndNodes, runwaySurface, runwayWidth, refs))

                # way: aeroway=taxiway
                elif tags['aeroway'] == 'taxiway':
                    print '\nFound a taxiway.  Tags:\n', tags

                    # Check to see if the taxiway is abandonded, if so skip it.
                    if isAbandoned(tags):
                        print 'NOTICE: Taxiway is abandoned, skipping. Way ID: ', osmid
                        continue

                    # TODO: probably need to do something better here, WED requires a name for all taxiway segments this default value of 'T' is just a hack for now to allow me to test them in game.
                    name = coalesceValue(('ref', 'name'), tags, 'T')

                    # Determine the taxiway surface type or assume 'concrete' if not specified.
                    if 'surface' in tags:
                        surface = tags['surface']
                    else:
                        surface = 'asphalt'

                    # Determine the taxiway width or use a sensible default.
                    if 'width' in tags:
                        width = tags['width']
                    else:
                        width = '15 m'

                    # We have successfully read all the data for this taxiway
                    # so add it to the list of completed taxiways to be put
                    # in the output file.
                    self.taxiways.append(Taxiway(name, surface, width, refs))

                # way: aeroway=apron
                if tags['aeroway'] == 'apron':
                    print '\nFound an apron.  Tags:\n', tags

                    name = coalesceValue(('name', 'ref'), tags, '')
                    surface = coalesceValue(('surface'), tags, 'concrete')

                    self.aprons.append(Apron(name, refs, surface))

    # Loops over all the object lists and tries to build a Shapely geometry
    # object for eact object in the list.
    def buildGeometries(self):
        for id, lon, lat in self.coords:
            self.coordDict[id] = (lon, lat)

        for id, tags, coord in self.nodes:
            self.nodeDict[id] = (tags, coord)

        # Loop over each list of objects stored in the Osm2apt class.
        listOfLists = self.objectLists
        for ls in listOfLists:
            for i in ls:
                i.buildGeometry(self.coordDict, self.nodeDict)

        for i in self.aerodromes:
            i.buildGeometry(self.coordDict, self.nodeDict)

    # Loops over all the objects in each object list (except for the aerodromes
    # list) and tries to place the given object with the closest areodrome it
    # can find (up to a maximum threshold distance).  If the object is
    # assosciated with an aerodrome, it is removed from the corresponding
    # Osm2apt list and placed in the relevant list for the assosciated
    # aerodrome.  Therefore, after this function finishes, the self.runways,
    # etc, lists should be empty, any objects they contain after that point
    # will not show up in X-plane as they will never be written to the apt.dat
    # file.
    def assosciateObjects(self):
        objectsToRemove = []

        if len(self.aerodromes) < 1:
            return

        for ls in self.objectLists:
            for obj in ls:
                shortestDistance = self.aerodromes[0].geometry.distance(obj.geometry)
                nearestAerodrome = self.aerodromes[0]
                # TODO: Refactor this and other occurrences into a function.
                for a in self.aerodromes:
                    tempDistance = a.geometry.distance(obj.geometry)
                    if tempDistance <= shortestDistance:
                        shortestDistance = tempDistance
                        nearestAerodrome = a

                # If the nearest aerodrome is within 0.02 degrees (about 2 km)
                # then assign the object to that aerodrome and remove it from
                # the general list.
                if shortestDistance < 0.02:
                    objectsToRemove.append((ls, obj))
                    nearestAerodrome.assosciatedObjects.append(obj)
                    obj.parentAerodrome = nearestAerodrome

        for ls, obj in objectsToRemove:
            ls.remove(obj)

        for a in self.aerodromes:
            a.cleanGeometries()

# Main function
print 'osm2apt version 0.4.0'

# Check and parse commandline
argumentParser = argparse.ArgumentParser()
argumentParser.add_argument('inputFile')
args = argumentParser.parse_args()
inputFileName = args.inputFile

print '=== Parsing OSM File ==='
osm2apt = Osm2apt_class()

parser = OSMParser(coords_callback=osm2apt.coordsCallback, nodes_callback=osm2apt.nodesCallback, ways_callback=osm2apt.waysCallback)
parser.parse(inputFileName)

print '\n\n\n=== Results ==='

print 'Successfully read in %s aerodromes\n' % len(osm2apt.aerodromes)
print 'Successfully read in %s runways' % len(osm2apt.runways)
print 'Successfully read in %s taxiways' % len(osm2apt.taxiways)
print 'Successfully read in %s aprons' % len(osm2apt.aprons)
print 'Successfully read in %s windsocks' % len(osm2apt.windsocks)
print 'Successfully read in %s beacons' % len(osm2apt.beacons)
print 'Successfully read in %s parkingPositions' % len(osm2apt.parkingPositions)
print 'Successfully read in %s lighted objects' % len(osm2apt.lightedObjects)

print '\nRead in %s osm nodes' % len(osm2apt.coords)

print '\n\n\n=== Building Geometries ==='
osm2apt.buildGeometries()
print 'Done.'

print '\n\n\n=== Assosciating Objects With Aerodromes ==='
osm2apt.assosciateObjects()
print 'Done\n'
print 'After trying to assosciate each read in object to a nearby aerodrome\n\
there are this many remaining un-assosciated objects (all these numbers should\n\
be 0, as any un-assosciated objects will not be put in the output apt.dat file):'
print "\t%s\tRunways" % len(osm2apt.runways)
print "\t%s\tTaxiways" % len(osm2apt.taxiways)
print "\t%s\tAprons" % len(osm2apt.aprons)
print "\t%s\tWindsocks" % len(osm2apt.windsocks)
print "\t%s\tBeacons" % len(osm2apt.beacons)
print "\t%s\tParking positions" % len(osm2apt.parkingPositions)
print "\t%s\tLighted Objects" % len(osm2apt.lightedObjects)

# TODO: Add all unassociated objects in these lists to the overpass query.

print '\n\n\n=== Outputing apt.dat ==='
outputFile = open('apt.dat', 'w')
outputFile.write('I\n');
outputFile.write('1000 Version - Imported data from OpenStreetMap processed by osm2apt.py - OpenStreetMap data is licenced under the terms of the ODBL\n');

for aerodrome in osm2apt.aerodromes:
    outputFile.write('\n\n')
    outputFile.write(aerodrome.toString())

outputFile.write('\n\n99\n\n')
overpassQueryFile.write(');\n\nout meta;\n')
print 'Done.'
print 'apt.dat written successfully.'

