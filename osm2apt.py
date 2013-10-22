#!/usr/bin/env python

import sys

from imposm.parser import OSMParser
from shapely.geometry import *

def lookahead(iterable):
    it = iter(iterable)
    last = it.next()
    for val in it:
        yield last, False
        last = val
    yield last, True

def metersToDeg(meters):
    return (meters / (6371000 * 2 * 3.1415927)) * 360

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

def printArea(area):
    ret = ''
    if area.exterior.is_ccw:
        coords = area.exterior.coords
    else:
        coords = reversed(area.exterior.coords)

    for coord, isLast in lookahead(coords):
        if not isLast:
            ret += '111 {0} {1}\n'.format(coord[1], coord[0])
        else:
            ret += '113 {0} {1}\n'.format(coord[1], coord[0])

    return ret

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

    def buildGeometry(self, coordDict):
        return

class Aerodrome(SpatialObject):

    def __init__(self, name, code, ele, nodes):
        self.name = name
        self.code = code
        self.ele = convertToUnit(ele, 'ft')
        self.nodes = nodes
        self.assosciatedObjects = []

    def buildGeometry(self, coordDict):
        if checkNodesFormClosedWay(self.nodes):
            self.geometry = Polygon(nodesToCoords(self.nodes, coordDict))
        else:
            print 'Not building geometry for aerodrome since it is not closed.'

    def toString(self):
        # Print out the main airport header line
        tempString = '1 {0} 0 0 {1} {2}\n'.format(self.ele, self.code, self.name)

        # TODO: Print out the boundary area of the aerodrome.

        # Loop over all of the taxiways and determine the unique set
        # (automatically has no duplicates since it is a python set) of nodes
        # used for the geometry of the taxiway network.
        # TODO: This should also include nodes for the runways when the on-runway taxiways are added.
        # TODO: If instead of just being a set this were a dictionary mapping nodes onto lists of ways this set could serve the same purpose but also be used to explor the topology of the taxiway network to find things like holding points, etc.
        taxiways = []
        taxiwayCoords = set()
        for obj in self.assosciatedObjects:
            tempString += obj.toString()
            if isinstance(obj, Taxiway):
                taxiways.append(obj)

        # If the airport has associated taxiways, also print out a taxiway network.
        if len(taxiways) > 0:
            # Start by building a set of the nodes used for taxiways at this specific airport.
            for taxiway in taxiways:
                for coord in taxiway.taxiwayCoords:
                    taxiwayCoords.add(coord)

            # Now print out the lines for all of the taxiway nodes for this airport.
            tempString += '1200\n'
            for coord in taxiwayCoords:
                tempString += '1201 {0} {1} both {2}\n'.format(coord[1], coord[0], coord[2])

            # Finally, print out the actual edges between the nodes.
            for taxiway in taxiways:
                prevCoord = 0
                for coord in taxiway.taxiwayCoords:
                    if prevCoord == 0:
                        prevCoord = coord
                    else:
                        # TODO: Implement oneway taxiways.
                        tempString += '1202 {0} {1} twoway taxiway {2}\n'.format(prevCoord[2], coord[2], taxiway.name)
                        prevCoord = coord
                    
        return tempString

class Runway(SpatialObject):

    def __init__(self, runwayEndNames, runwayEndNodes, surface, width, nodes):
        self.runwayEndNames = runwayEndNames
        self.runwayEndNodes = runwayEndNodes
        self.surface = surface
        self.surfaceInteger = surfaceStringToInt(self.surface)
        self.width = convertToUnit(width, 'm');
        self.nodes = nodes

    def buildGeometry(self, coordDict):
        self.geometry = LineString(nodesToCoords(self.nodes, coordDict))
        # Check to see if the first end of the runway is drawn with the lower
        # numbered end first or if it needs to be reversed so that the
        # direction of the line lines up with the name.  Note that this assumes
        # the ref tag gives the lower number first, like 18/36.
        #TODO: Need to fix this runway reversing code, for now it seems to not be correct so it is commented out.
        '''
        if ((self.geometry.coords[0][0] < self.geometry.coords[-1][0]) or ((self.geometry.coords[0][0] == self.geometry.coords[-1][0]) and self.geometry.coords[0][1] > self.geometry[-1][1])):
            self.geometry.coords = list(self.geometry.coords)[::-1]
            self.nodes.reverse()
            print 'Reversing runway %s' % self.runwayEndNames
            print list(self.geometry.coords)
        else:
            print 'Not reversing runway %s' % self.runwayEndNames
            print list(self.geometry.coords)
        '''

    def toString(self):
        return '100 {0} {1} 0 0.15 0 0 1 {2} {3} {4} 0 0 1 0 0 0 {5} {6} {7} 0 0 1 0 0 0\n'.format(self.width, self.surfaceInteger, self.runwayEndNames[0], self.geometry.coords[0][1], self.geometry.coords[0][0], self.runwayEndNames[1], self.geometry.coords[-1][1], self.geometry.coords[-1][0])

class Taxiway(SpatialObject):

    def __init__(self, name, surface, width, nodes):
        self.name = name
        self.surface=surface
        self.surfaceInteger = surfaceStringToInt(self.surface)
        self.width = convertToUnit(width, 'm');
        self.nodes = nodes

    def buildGeometry(self, coordDict):
        self.taxiwayCoords = nodesToCoords(self.nodes, coordDict)
        self.geometry = LineString(self.taxiwayCoords)
        self.concreteGeometry = self.geometry.buffer(metersToDeg(self.width), 2)

    def toString(self):
        ret = '110 {0} 0.15 360 {1}\n'.format(self.surfaceInteger, self.name)
        ret += printArea(self.concreteGeometry)

        return ret

class Windsock(SpatialObject):

    def __init__(self, coord, lit):
        self.coord = coord
        if lit in ('yes', 'true', '1'):
            self.lit = 1
        else:
            self.lit = 0

    def buildGeometry(self, coordDict):
        self.geometry = Point(self.coord)

    def toString(self):
        return '19 {0} {1} {2} WS\n'.format(self.coord[1], self.coord[0], self.lit)

class Apron(SpatialObject):

    def __init__(self, name, nodes, surface):
        self.name = name
        self.nodes = nodes
        self.surface = surface
        self.surfaceInteger = surfaceStringToInt(self.surface)

    def buildGeometry(self, coordDict):
        if checkNodesFormClosedWay(self.nodes):
            # TODO: Need to make sure the coords in the resulting geometry for a counter clockwise ring.
            self.geometry = Polygon(nodesToCoords(self.nodes, coordDict))
        else:
            print 'Not building geometry for apron since it is not closed.'

    def toString(self):
        ret = '110 {0} 0.15 360 {1}\n'.format(self.surfaceInteger, self.name)
        ret += printArea(self.geometry)

        # If this apron is named, also create a start location for it.
        # TODO: If an airport has no named aprons we should choose the largest one and put a start location there.  An alternative is for the aerodrome class to make sure that it has at least one named apron by naming the largest on 'Ramp' or something if none are named.
        if len(self.name) > 0:
            startupLoc = self.geometry.representative_point()
            # TODO: Set the startup heading to be something reasonable, toward the nearest taxiway would be a good one.
            ret += '15 {0} {1} 360 {2}\n'.format(startupLoc.y, startupLoc.x, self.name)

        return ret
class Osm2apt_class(object):

    aerodromes = []

    # Every one of the lists in this next section should be added to
    # objectLists like 'runways' is.
    objectLists = []
    runways = [];    objectLists.append(runways)
    taxiways = [];   objectLists.append(taxiways)
    windsocks = [];   objectLists.append(windsocks)
    aprons = [];   objectLists.append(aprons)

    coords = []
    coordDict = {}

    # Callback method to simply read in all node ID's and coordinates
    def coordsCallback(self, coords):
        for c in coords:
            self.coords.append(c)

    # Callback method to read all nodes which have tags on them
    def nodesCallback(self, nodes):
        for osmid, tags, coord in nodes:
            if 'aeroway' in tags:
                #node: aeroway=windsock
                if tags['aeroway'] == 'windsock':
                    print 'Windsock'
                    lit = coalesceValue(('lit'), tags, 'no')
                    self.windsocks.append(Windsock(coord, lit))

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
                    if 'icao' in tags:
                        aerodromeCode = tags['icao']
                    elif 'faa' in tags:
                        aerodromeCode = tags['faa']
                    else:
                        print 'ERROR: Aerodrome does not have an ICAO or FAA code, skipping.  Way ID: ', osmid
                        continue

                    print 'Airport Code is: ', aerodromeCode

                    # Get the airport elevation, if there is no 'ele' tag skip
                    # the airport.
                    if 'ele' in tags:
                        aerodromeELE = tags['ele']
                    else:
                        print 'ERROR: Aerodrome does not have an elevation (ele), skipping. Way ID: ', osmid
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
                    if 'ref' in tags:
                        runwayEndNames = tags['ref'].split('/')
                        print 'Runway end names are: ', runwayEndNames
                        print 'Runway end nodes are: ', runwayEndNodes
                    else:
                        print 'ERROR: Runway does not have a "ref" tag set, skipping.  Way ID: ', osmid
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

                    #TODO: Handle displaced throsholds and blast pads.

                    #TODO: Handle runway markings, TDZ lighting, and REIL lighting.

                    # We have successfully read all the data for this runway
                    # so add it to the list of completed runways to be put
                    # in the output file.
                    self.runways.append(Runway(runwayEndNames, runwayEndNodes, runwaySurface, runwayWidth, refs))

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

        # Loop over each list of objects stored in the Osm2apt class.
        listOfLists = self.objectLists
        for ls in listOfLists:
            for i in ls:
                i.buildGeometry(self.coordDict)

        for i in self.aerodromes:
            i.buildGeometry(self.coordDict)

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
                for a in self.aerodromes:
                    tempDistance = a.geometry.distance(obj.geometry)
                    if tempDistance <= shortestDistance:
                        shortestDistance = tempDistance
                        nearestAerodrome = a

                # If the nearest aerodrome is within 0.05 degrees (about 5 km)
                # then assign the object to that aerodrome and remove it from
                # the general list.
                if shortestDistance < 0.05:
                    objectsToRemove.append((ls, obj))
                    nearestAerodrome.assosciatedObjects.append(obj)

        for ls, obj in objectsToRemove:
            ls.remove(obj)

# Main function
print 'osm2apt version 0.1.0'

'''
print '=== Testing conversion ratios ==='
testStrings = ['7', '7   ', '   7', '7 ft', ' 7 ft', '7     ft', '7 ft   ', '7 m', '3 blobs']
for string in testStrings:
    print '"%s" is %f m' % (string, convertToUnit(string, 'm'))
print convertToUnit('7 ft', 'hello')
print '\n\n\n'
'''

# Check and parse commandline
#TODO: Imlement this.

print '=== Parsing OSM File ==='
osm2apt = Osm2apt_class()

parser = OSMParser(coords_callback=osm2apt.coordsCallback, nodes_callback=osm2apt.nodesCallback, ways_callback=osm2apt.waysCallback)
parser.parse('Airports.osm')

print '\n\n\n=== Results ==='

print 'Successfully read in %s aerodromes\n' % len(osm2apt.aerodromes)
print 'Successfully read in %s runways' % len(osm2apt.runways)
print 'Successfully read in %s taxiways' % len(osm2apt.taxiways)
print 'Successfully read in %s aprons' % len(osm2apt.aprons)
print 'Successfully read in %s windsocks' % len(osm2apt.windsocks)

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

print '\n\n\n=== Outputing apt.dat ==='
outputFile = open('apt.dat', 'w')
outputFile.write('I\n');
outputFile.write('1000 Version - Imported data from OpenStreetMap processed by osm2apt.py - OpenStreetMap data is licenced under the terms of the ODBL\n');

for aerodrome in osm2apt.aerodromes:
    outputFile.write('\n\n')
    outputFile.write(aerodrome.toString())

outputFile.write('\n\n99\n\n')
print 'Done.'
print 'apt.dat written successfully.'

