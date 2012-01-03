#!/usr/bin/env python

##     GPS-data clustering software
##     Copyright (C) 2008 mike warren (mike@mike-warren.com)
##
##     This program is free software: you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation, either version 3 of the License, or
##     (at your option) any later version.
##
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with this program.  If not, see <http://www.gnu.org/licenses/>.

##
## this takes as input a CSV file of GPS-collar fixes. A typical such
## file looks like this:
##
## FIX_#,CASE,GMT_DATE,GMT_TIME,LMT_DATE,LMT_TIME,N7,N8,N9,LATITUDE,LONGITUDE,HEIGHT,N13,N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,N36,N37,N38,N39,N40,N41,N42,N43,N44,N45,N46
## 2,1,2/18/2006,22:04:18,2/18/2006,22:04:18,-1710094,-3498899,5035956,52.47025,-116.04717,1396.71,6.4,3D,Yes,5,15,41,16,41,21,35,26,38,29,35,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.42,3.56,11,N/A,N/A,N/A
## 4,1,2/19/2006,4:01:24,2/19/2006,4:01:24,-1710121,-3498898,5035938,52.47007,-116.04753,1389.11,4,3D,Yes,5,1,41,11,41,14,38,20,38,25,35,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.42,3.62,19,N/A,N/A,N/A
##
## The file-reading code presumes the first line of the file is
## column-headers and discards it.
##
## If you wish to change how the data is parsed, see the constructor
## for FixPoint, which gets a tuple containing all the data from one
## line. Currently, it only parses out the following columns:
##
##    0: fix number
##    1: "case"
##    4 + 5: time of the fix
##    9 + 10: latitude and longitude
##    11: height (actually not used currently)
##
## There are a couple of options you may want to change; see "OPTIONS"
## below.
##



import traceback
import csv
import math
import time
import sys
import os
import os.path


##
## OPTIONS you might want to change
##

# approximate target cluster radius, in meters (depends a little on
# the method selected below).
CLUSTER_RADIUS = 25
#ATF changed this from 200 to 25 following augstine 2003, "Spatial Heterogeneity in the Herbaceous Layer of a Semi-Arid Savanna Ecosystem"
#ATF changed CLUSTERMETHODS to 'centroid'
# available clustering methods; see around line 200 for more
# explanation
CLUSTER_METHODS = ['centroid']
CLUSTER_METHOD = CLUSTER_METHODS[0]

# debugging stuff
GOOGLE_MAP = False                      # output google-map data (to plot on a google-map)
DEBUG_DISTANCE_MATRIX = False           # output a matrix of all points with the distance to all others
EXPORT_RAW = False                      # pickle the clusters to "clusters.pickle"
LOAD_PICKLED_CLUSTERS = False           # don't re-cluster; just load whatever is in "clusters.pickle"


## code follows





#hash table of FixPoint objects indexed by fix number
dataByFixNumber = {}

def greatCircleDistance(a,b):
    """ from http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/393241
    much more info at: http://en.wikipedia.org/wiki/Great_circle_distance

    this is using the haversine arcsin method (which won't work
    well for points on opposite sides of the earth, but that
    shouldn't be a problem for this data)

    this will return the distance in METERS

    "other" can be another FixPoint or a 2-tuple of (lat,lng)

    """


    (alat,alng) = a
    (blat,blng) = b

    lngdist = math.radians(blng - alng)
    latdist = math.radians(blat - alat)

    alat = math.radians(alat)
    blat = math.radians(blat)
    alng = math.radians(alng)
    blng = math.radians(blng)
    
    a = (math.sin(latdist / 2))**2 + math.cos(alat) * math.cos(blat) * (math.sin(lngdist / 2))**2
    c = 2 * math.asin(min(1, math.sqrt(a)))
    # dist = 6399.592 * c  # magic number for poles
    # dist = 6335.437 * c  # magic number for equator
    dist = 6372.795 * c # average arcradius of Earth (in km)
    return int(dist * 1000)              # return meters, not kilometers (don't return decimals; up to 0.5% error anyway)


##
## each of the points in the file are read as one of these
## if there was no data for that point, self.latlng will be None
##

class FixPoint:
    def __init__(self, *args):
        self.number = int(args[0])
        self.case = args[1]
        
        # excel seems to do this:
        self.time = time.mktime(time.strptime(args[4] + ' ' + args[5], "%m/%d/%Y %H:%M:%S"))
        
        # gnumeric seems to do this:
        if 0:
            thetm = args[4] + ' ' + args[5]
            if thetm[-7:] == '.000000':
                thetm = thetm[:-7]
            self.time = time.mktime(time.strptime(thetm, "%Y/%m/%d %H:%M:%S"))

        self.lmt = args[4] + ' ' + args[5]
        
        if args[9] == 'N/A' or args[10] == 'N/A':
            self.latlng = None
        else:
            self.latlng = (float(args[9]), float(args[10]))
        #self.height = float(args[11])

    def distance(self,other):
        """distance from me to the other point, either a (lat,lng)
        tuple or another FixPoint"""
        
        if self.latlng is None:
            return None
        if isinstance(other,FixPoint) and other.latlng is None:
            return None

        if not isinstance(other,FixPoint):
            (olat,olng) = other
        else:
            (olat,olng) = other.latlng
            

        return greatCircleDistance( self.latlng, (olat,olng) )


    def timedelta(self,other):
        """
        return the (absolute) difference in time between this and the other point in HOURS
        """

        if 0:
            # return time difference based on fix number (3 hours apart), atf changed to 600, as in 600 secs = 10min
            return math.fabs( (self.number*600) - (other.number*600) )

        else:
            # actual seconds-based time difference
            seconds = other.time - self.time
            seconds = math.fabs(seconds)
            return seconds            # return hours ATF changed to nothing, returns seconds.

    def __str__(self):
        if self.latlng:
            return "<FixPoint #%d (%f,%f) %f seconds>" % (self.number,self.latlng[0],self.latlng[1],self.time)
        return "<FixPoint #%d (N/A) %f seconds>" % (self.number,self.time)

    def __repr__(self):
        return self.__str__()

##
## This represents a cluster, obviously
##

#########ATF CHANGED BELOW: changed (self.points[-1].timedelta(a) <= 24*6) or (self.points[0].timedelta(a) <= 24*6):
###if (self.points[-1].timedelta(a) <= 3600) or (self.points[0].timedelta(a) <= 3600):
   ##             return True
     ##       if (self.points[-1].timedelta(b) <= 1.0) or (self.points[0].timedelta(b) <= 1.0):
       ##         return True    

class Cluster:
    def __init__(self):
        self.points = []
        self.centroidcache = None

    def clusterInCluster(self,c):
        """returns True if the other cluster c is"inside" this one
        (i.e. if they can be merged); centroids must be within CLUSTER_RADIUS
        and the points must be within 6 days
        """
        
        if greatCircleDistance(self.centroid(),c.centroid()) < CLUSTER_RADIUS:
            a = c.points[0]
            b = c.points[-1]
            if (self.points[-1].timedelta(a) <= 3600) or (self.points[0].timedelta(a) <= 3600):
                return True
            if (self.points[-1].timedelta(b) <= 3600) or (self.points[0].timedelta(b) <= 3600):
                return True
        return False
    
    def allFixNumbers(self):
        rtn = ''
        for x in self.points:
            rtn = rtn + '%d ' % x.number
        return rtn

    def pointInCluster(self,p):
        if p.latlng is None:
            return False

###ATF changed below "<= 6*24.0" to "<= 1.0"
        
        if len(self.points) == 0:
            return False
        elif len(self.points) == 1:
            timeok = (self.points[-1].timedelta(p) <= 3600) or (self.points[0].timedelta(p) <= 3600)
            return (self.points[0].distance(p) < CLUSTER_RADIUS) and timeok
        else:
            timeok = (self.points[-1].timedelta(p) <= 3600) or (self.points[0].timedelta(p) <= 3600)
            if CLUSTER_METHOD == 'any':
                # method "is new point within CLUSTER_RADIUS m of ANY point currently in the cluster"?
                for x in self.points:
                    if p.distance(x) <= CLUSTER_RADIUS:
                        if timeok:
                            return True
                return False

            if CLUSTER_METHOD == 'all':
                # method "is new point within CLUSTER_RADIUS m of ALL points currently in the cluster"?
                for x in self.points:
                    if p.distance(x) > CLUSTER_RADIUS:
                        if timeok:
                            return False
                return True

            if CLUSTER_METHOD == 'centroid':
                # method "is new point within CLUSTER_RADIUS m of current centroid of cluster"?
                return p.distance(self.centroid()) < CLUSTER_RADIUS #and self.points[-1].timedelta(p) < 1.0 <- ATF CHANGED COMMENT to 1.0 from 24.0

        # default, not in cluster
        return False
        
    def addPoint(self,p):
        if p.latlng is None:
            print "ERROR",p
            return
        if p not in self.points:
            self.points.append(p)
        self.points.sort( lambda x,y: cmp(x.number,y.number) )

    def centroid(self):
        if self.centroidcache is not None and self.centroidcache[1] == len(self.points):
            return self.centroidcache[0]
        avglat = 0.0
        avglng = 0.0
        for p in self.points:
            (lat,lng) = p.latlng
            avglat = avglat + lat
            avglng = avglng + lng
        rtn = (avglat/len(self.points), avglng/len(self.points))
        self.centroidcache = (rtn,len(self.points))
        return rtn

    def averageDistanceFromCentroid(self):
        avgdist = 0.0
        for x in self.points:
            avgdist = avgdist + x.distance(self.centroid())
        return avgdist / len(self.points)

    def maxDistanceFromCentroid(self):
        max = 0.0
        for x in self.points:
            dist = x.distance(self.centroid())
            if dist > max:
                max = dist
        return max
        

    def totalTime(self):
        """time spread of cluster in hours"""
        for x in range(1,len(self.points)):
            assert self.points[x].number > self.points[x-1].number
        return self.points[-1].timedelta(self.points[0])

    def totalTheoreticalPoints(self):
        """
        return the number of points which "should" be in the cluster; this is the difference between the biggest and
        smallest fix-number in the cluster
        """
        
        assert len(self.points) > 1
        return self.points[-1].number - self.points[0].number + 1

    def missingPoints(self):
        """how many points are "missing" from the cluster?
        that is, if the first point is #10 and the last is
        #20, there "should" be 11 points -- but if, say,
        #12, #13 and #14 aren't in this cluster [and they might be missing
        totally -- i.e. "N/A" fixes] but the rest are, then this method should return 3.
        """
        
        last = self.points[0]
        missing = 0
        for x in self.points[1:]:
            if last.number != x.number-1:
                assert (x.number-last.number-1) > 0
                missing = missing + (x.number-last.number-1)
            last = x
        return missing

    def fixesAwayFromCluster(self,returnListOfMissing=False):
        """
        similar to above, but say that #13 in the
        above example was a "successful" fix (i.e.
        not "N/A" in the file) then this method
        should return 1.
        """

        last = self.points[0]
        away = 0
        bonus = []
        for x in self.points[1:]:
            if last.number != x.number-1:
                for n in range(last.number+1, x.number):
                    if dataByFixNumber.has_key(n):
                        assert(n not in self.points)
                        fixn = dataByFixNumber[n]
                        if fixn.latlng != None:
                            away = away + 1
                            bonus.append(n)
            last = x
        if returnListOfMissing:
            return (away,bonus)
        return away

    def numberOfhours(self):
        """
        Kyle: The number of 24 hour periods where a fix was obtained
        at the cluster.  For example, a cluster might have 8 points
        and span 4 days, but 7 of the points were obtained on the
        first day and one was obtained on the last day.  This cluster
        has two 24 hour periods where there was a fix. A similar 8
        point cluster has 2 points obtained on each of 4 days and so
        has a score of 4 in this output.

        Mike: Algorithm: starting with the first point, we "quantize"
        everything into 24-hour-long buckets, starting with the hour
        of the first fix (i.e. if there was a fix at 11:30pm and two
        more at 2am, they'd only be on one"day")
        """
###ATF changed 24*60*60 to (60*60) so is by hour (3600 secs) not day buckets....
        hourstarts = []
        for x in self.points:
            if len(hourstarts) == 0 or x.time > hourstarts[-1] + (60*60):
                hourstarts.append(x.time)
        return len(hourstarts)
            

    def __eq__(self,other):
        """return True if this cluster is equal to "other" -- has all the same points"""
        
        if len(self.points) == 0 or len(other.points) == 0:
            return False
        if len(self.points) != len(other.points):
            return False
        for x in self.points:
            if not x in other.points:
                return False
        return True

    def __str__(self):
        if True:
            pnt = '['
            for x in self.points:
                pnt = pnt + str(x.number) + ', '
            pnt = pnt[:-2] + ']'
            return 'Cluster of %d points (%d hours spread, missing %d points, %dm average distance).' % (len(self.points), int(self.totalTime()), self.missingPoints(),int(self.averageDistanceFromCentroid())) + pnt

        else:
            rtn = 'Cluster of %d points (%f hours spread, missing %d points, %dm average distance):\n' % (len(self.points), self.totalTime(), self.missingPoints(),int(avgdist))
            start = self.points[0]
            avgdist = self.averageDistanceFromCentroid()
            for x in self.points:
                rtn = rtn +  "   #%d (%f,%f) timedelta %f hours\n" % (x.number,x.latlng[0],x.latlng[1],x.timedelta(start))
            rtn = rtn + '   average distance from centroid: %dm\n' % int(avgdist)
            return rtn

    def html(self):
        rtn = 'Cluster of %d points<br/>%f hours spread.<br/>missing %d points<br/>' % (len(self.points), self.totalTime(), self.missingPoints())
        rtn = rtn + 'avg dist from centroid: %dm' % int(self.averageDistanceFromCentroid())
        return rtn




def processFile(fname):
    print 'Processing "%s"' % fname
    
    # open file, read all lines and trim off the first one (column headers)
    try:
        lines = open(fname).readlines()[1:]
    except IOError:
        print 'Error opening/reading "%s"...' % fname
        return

    # parse it as a CSV file
    input = csv.reader(lines)

    # "data" is any data with actual coordinates
    # "dead" is "N/A" points
    data = []
    dead = []

    try:
        last = None
        while True:
            foo = input.next()
            pt = FixPoint(*foo)
            dist = 0
            if last:
                dist = last.distance(pt)
                if dist is None: dist = -1
            if 1:
                print "  loaded fix #%d (%dm away)" % (pt.number,dist)
            last = pt
            if pt.latlng:
                data.append(pt)
            else:
                dead.append(pt)
            dataByFixNumber[pt.number] = pt
    except StopIteration:
        pass

    print "  loaded %d points: %d good, %d missing." % (len(data)+len(dead),len(data),len(dead))


    if True:
        print '  hierarchical clustering.'
        ## do "hierarchical clusering"
        ##
        ## this works by making a "cluster" for each point containing just that point.
        ## then, an attempt is made to merge each cluster with all other clusters until
        ## no more merging can be done.
        clusters = []

        ## data now pre-filtered to only have valid fixes
        
        ## make all points into a "cluster" of one point
        for x in data:
            c = Cluster()
            c.addPoint(x)
            clusters.append(c)


        ## merge clusters which match distance criterion (e.g. centers
        ## are within CLUSTER_RADIUS m of each other)

        print "  merging clusters.",
        sys.stdout.flush()
        def merge(clusters):
            for c in clusters:
                for n in clusters:
                    if c is n:
                        continue
                    if c.clusterInCluster(n):
                        print "   making new cluster:"
                        print "  ",c
                        print "  ",n
                        print "-----------"
                        foo = Cluster()
                        bar = Cluster()
                        for y in c.points:
                            foo.addPoint(y)
                        for y in n.points:
                            if foo.pointInCluster(y):
                                foo.addPoint(y)
                            else:
                                bar.addPoint(y)
                                
                        if ( (c == bar) and (n == foo)) or (( (c == foo) and (n == bar))):
                            print "new clusters same as old; not replacing"
                            continue
                        
                        clusters.remove(c)
                        clusters.remove(n)
                        clusters.append(foo)
                        if len(bar.points) > 0:
                            print "   2 new clusters:"
                            print foo
                            print "-----------"
                            print bar
                            print "-----------"
                            clusters.append(bar)
                        else:
                            print "   new cluster:"
                            print foo
                            print "-----------"
                        return True
            return False

        if not LOAD_PICKLED_CLUSTERS:
            last = time.time()
            while merge(clusters):
                now = time.time()
                if now - last > 2:
                    #print '    \b\b\b\b\b%04d' % len(clusters),
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    last = now
            print 'done\n  filtering:'

        else:
            import pickle
            clusters = pickle.load(open('clusters.pickle','r'))

        ## from the list of all clusters, filter out the
        ## ones we don't want

                                            # only keep clusters with > 1 point
        before = len(clusters)
        clusters = filter(lambda x: len(x.points) > 1, clusters)
        print "    %d 1-point clusters removed" % (before-len(clusters),)

        print "  found %d clusters:" % len(clusters)

        if EXPORT_RAW:
            import pickle
            pickle.dump(clusters,open('clusters.pickle','w'))

        if 0:
            ## debugging; if we can still merge here, something is wrong
            print "  more merging?",merge(clusters)

        if True:
            ## write cluster data to a .csv file
            csvfile = open(os.path.splitext(fname)[0] + '-mikeclusters.csv','w')
            csvfile.write('''"first fix number","first fix LMT","last fix number","theoretical total points","actual number of points","points away","fix success","time span (hours)","hour-periods","average distance (m)","cluster radius (m)","centroid latitude","centroid longitude","all fix numbers in cluster","all away points"\n''')
            for x in clusters:
                totalpoints = x.totalTheoreticalPoints()
                assert (x.missingPoints() >= 0)
                awayfixes, awaypoints = x.fixesAwayFromCluster(True)
                # make awaypoints into a space-separated string for output
                awaypoints = ' '.join(map(str,awaypoints))
                
                fixsuccess = ((float(totalpoints) - (totalpoints - (awayfixes+len(x.points)))) / float(totalpoints)) * 100.0
                hours = x.numberOfhours()
                csvfile.write('%d,%s,%d,%d,%d,%d,%f,%d,%d,%d,%d,%f,%f,%s,%s\n' % (x.points[0].number,x.points[0].lmt,x.points[-1].number,totalpoints,len(x.points),awayfixes,fixsuccess,x.totalTime(),hours,int(x.averageDistanceFromCentroid()),int(x.maxDistanceFromCentroid()),x.centroid()[0],x.centroid()[1],x.allFixNumbers(),awaypoints))
            csvfile.close()

        if GOOGLE_MAP:
            ## write data to produce the google-map
            file = open('coords','w')
            last = None
            for p in data:
                if p.latlng:
                    if last and last.number != p.number-1:
                        file.write('geo:lat=%f geo:lon=%f end\n' % p.latlng)
                    else:
                        file.write('geo:lat=%f geo:lon=%f waypoint\n' % p.latlng)
                last = p

            for x in clusters:
                print "  ",x
                lat,lng = x.centroid()
                file.write('geo:lat=%f geo:lon=%f %s\n' % (lat,lng,x.html(),))

            file.close()

        return clusters
        

            


    if DEBUG_DISTANCE_MATRIX:
        ## this outputs a NxN matrix of each point and the distance to all other points
        ## (for debugging)
        print 'writing "matrix.csv"'
        f = open('matrix.csv','w')
        for x in data:
            for y in data:
                dist = x.distance(y)
                if dist is None: dist = -1
                f.write('%d,' % dist)
            f.write('\n')
        f.close()




if __name__ == "__main__":
    ## process all files from the command-line
    if len(sys.argv) > 1:
        for x in sys.argv[1:]:
            if 'mikeclusters' in x:
                continue
            processFile(x)

    else:
        for x in os.listdir('.'):
            if os.path.splitext(x)[1] == '.csv':
                if not 'mikecluster' in x:
                    try:
                        processFile(x)
                    except:
                        print "ERROR processing",x
                        for y in traceback.format_tb(sys.exc_info()[2]):
                            print y,
                        print "CONTINUING"
                else:
                    print x,"looks like a cluster-file; skipping"

            else:
                print 'Skipping "%s" because it doesn\'t end in .csv' % x

        raw_input('press return')
