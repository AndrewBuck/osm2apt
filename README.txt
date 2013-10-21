NOTE: This README includes much future planned functionality which is not
yet implemented in the actual script, it is designed to give an overview of the
design philosophy of the program and where it is going.  See the version info
at the end of this README for details on what has been implemented so far.

=== Introduction ===

osm2apt.py - A python script to convert an OSM .xml or .pbf file containing one
or more airfields and associated runways, taxiways, windsocks, etc, into the
Airport.dat file format used by the X-plane flight simulator (specifically
targeting X-plane 10, which is the current version at this time).  The OSM file
can be downloaded from Overpass, XAPI, or even the regular API and requires no
preprocessing or filtering (any data that is not specifically airport related
will be either ignored or used to help inform guesses about airport flows, etc,
by doing things like having a traffic pattern be over empty area at the edge of
town rather than over a residential area, etc).

The program reads in the OSM data, stores it into some python objects to allow
it to do internal manipulations of the data more easily than in the native OSM
format and finally writes out the resulting airport in the format that X-plane
expects.

The hope is to have these generated airports convey as much real world
information that exists in OSM as possible, while at the same time making
sensible guesses for any information not explicitely entered into OSM (such as
runway widths, taxiway surfaces, etc), and finally to output an airport that
takes full advantage of the features that X-plane 10 affords.  For example,
exported airports should have proper concrete surfaces for each taxiway, as
well as an assosciated taxiway network which is as routable as the original OSM
data, and all taxiways should be properly named if they are in OSM, etc; it is
also planned to make approximate airport flows and other "higher level" kinds
of airport information as best as possible based on the number and properties
of the runways and possibly other nearby OSM data to inform the flow
information.  After the apt.dat file is written it should either be useable in
X-plane directly or be improved with minimal tweaking in WED to touch up
cosmetic issues, but the bulk of the taxiway network, airport flows, surfaces,
and exclusion areas should already be made to a satisfactory level so that the
user does not have to spend a lot of time doing the 'grunt work' of making all
of these redundant things manually, which is both time consuming and error
prone if done manually.



=== Dependencies ===

*  The ImpOSM Parser (imposm.parser) - Included in the Git repo as a submodule.
*  Shapely
*  Numpy

git submodule update --init --recursive
sudo apt-get install python-pip
sudo pip install Shapely

=== Version Info ===

0.1.0 - Implementing very basic routines to just read a few object types
(aerodrome, runway, etc) and output a simple .dat representation which handles
just basic 'private' fields with very simple layouts.  Basically a 'proof of
concept' version.  Nothing done yet with taxiway networks, airport flows, or
exclusion areas.

