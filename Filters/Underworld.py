# This filter is for generating underworld features (Caves, tunnels, bottomless pits to the void, etc).
# Suggested by a number of people, most recently @lemoesh
# abrightmoore@yahoo.com.au
# http://brightmoore.net
# My filters may include code and inspiration from PYMCLEVEL/MCEDIT mentors @Texelelf, @Sethbling, @CodeWarrior0, @Podshot_

from pymclevel import TAG_String # @Texelelf
from pymclevel import TileEntity # @Texelelf

import time # for timing
from math import sqrt, tan, sin, cos, pi, ceil, floor, acos, atan, asin, degrees, radians, log, atan2
from random import *
from numpy import *
from pymclevel import alphaMaterials, MCSchematic, MCLevel, BoundingBox, MCMaterials
from mcplatform import *
from os import listdir
from os.path import isfile, join
import glob
from pymclevel import TAG_List, TAG_Byte, TAG_Int, TAG_Compound, TAG_Short, TAG_Float, TAG_Double, TAG_String, TAG_Long
# import from @Texelelf
from copy import deepcopy
import bresenham # @Codewarrior0
from random import Random # @Codewarrior0
import inspect # @Texelelf
from PIL import Image
import png

# GLOBAL
CHUNKSIZE = 16
MAXANGLES = 360
a = 2*pi/MAXANGLES
AIR = (0,0)
LAVA = (11,0)
coloursList = "15 7 8 0 6 2 10 11 3 9 13 5 4 1 14 12".split() # I like this sequence
#coloursList = "15 15 15 15 15 11 11 11 0".split()
colours = map(int, coloursList)
palette = []
palette.append( [ (1,0),(1,1), (1,2), (1,3), (1,4), (1,5), (1,6), (80,0) ] )
palette.append( [ (95,3),(22,0), (24,2), (19,0), (18,0), (3,0), (1,0), (80,0) ]	)
RAND2 = Random(42)
	
# Choose a colour pallette
print 'Pallette Selection'
gap = RAND2.randint(1,4)
for i in range(32):
	newCols = []
	baseIndex = RAND2.randint(0,len(colours))

	for j in range(8):
		newCols.append( (159, (int)(colours[(baseIndex+j*gap)%len(colours)]) ) )
	palette.append(newCols)
print 'Pallette Selection complete'


#COUNTTUNNELINVOKES = 0 # Global

# Filter pseudocode:
#

inputs = (
		("UNDERWORLD", "label"),
		("Operation:", ( "Tunnels",
						"Tunnels SMALL",
						"Tunnels BIG",
						"Tunnels: Brownian (Random)",
						"Tunnels: Brownian (Random) SMALL",
						"Tunnels: Brownian (Random) BIG",
						"Tunnel: Smooth Path",
						"Ores: User selection only",
						"Ores: Stone Variants - Plain",
						"Ores: Ore Spheres - Plain",
						"Ores: Ore Spheres - No Layers",
						"Ores: Ore layers - Stained Clay",
						"Ores: Twisting Seam",
						"Columns",
						"Volcano",
						"Limestone Cave",
						"Interference",
						"Interference: Variable Frequency",
						"Perlin 3D",
						"Perlin 2D",
						"Denseland",
						"Web",
						"Maze: Stone TBD"
						)),
		("Material:", alphaMaterials.Glowstone),
		("Block Type to Replace:", alphaMaterials.Stone),
		("Chance:",100),
		("Seed:", 0),
		("Cache:",True),
		("abrightmoore@yahoo.com.au", "label"),
		("http://brightmoore.net", "label"),
)

def underworld(level, box, options):
	method = "Underworld"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start

	COUNTTUNNELINVOKES = 0
	SEED = options["Seed:"]
	OPER = options["Operation:"]
	material = (options["Material:"].ID,options["Material:"].blockData)
	chance = options["Chance:"]
	MAXITER = 1

	if SEED == 0:
		SEED = randint(0,999999999999)
		MAXITER = randint(1,16)
		print MAXITER
	R = Random(SEED)
	print SEED
	(px,py,pz) = (centreWidth,centreHeight,centreDepth)
	#for i in xrange(0,MAXITER):
	
	if OPER == "Tunnels":
		tunnelTwisting(level, box, options, R,(px,py,pz),0,-1,-1,chance)
	elif OPER == "Tunnels SMALL":
		tunnelTwisting(level, box, options, R,(px,py,pz),500,-1,-1,chance)
	elif OPER == "Tunnels BIG":
		tunnelTwisting(level, box, options, R,(px,py,pz),10000,-1,-1,chance)
	elif OPER == "Tunnels: Brownian (Random) BIG":
		tunnelBrownian(level,box,options,R,(px,py,pz),5000)
	elif OPER == "Tunnels: Brownian (Random)":
		tunnelBrownian(level,box,options,R,(px,py,pz),3000)
	elif OPER == "Tunnels: Brownian (Random) SMALL":
		tunnelBrownian(level,box,options,R,(px,py,pz),1000)
	elif OPER == "Ores: User selection only":
		repopulateSphere(level,box,options,R,material,3,1,3,box.miny,box.maxy,chance) # User selection
	elif OPER == "Ores: Stone Variants - Plain":
		# repopulateSphere(level,box,options,R,material,3,1,3,box.miny,box.maxy,chance) # User selection
		repopulateSphere(level,box,options,R,(1,1),8,1,3,box.miny,box.maxy,chance) # Granite
		repopulateSphere(level,box,options,R,(1,3),8,1,3,box.miny,box.maxy,chance) # Diorite
		repopulateSphere(level,box,options,R,(1,5),8,1,3,box.miny,box.maxy,chance) # Andesite
		repopulateSphere(level,box,options,R,(3,0),8,1,3,box.miny,box.maxy,chance) # dirt
		repopulateSphere(level,box,options,R,(13,0),8,0,3,box.miny,box.maxy,chance) # Gravel
	elif OPER == "Ores: Ore Spheres - Plain":
		repopulateSphere(level,box,options,R,(14,0),2,1,1,1,29,chance) # Gold ore
		repopulateSphere(level,box,options,R,(15,0),2,1,1,1,64,chance) # Iron ore
		repopulateSphere(level,box,options,R,(16,0),3,1,1,box.miny,box.maxy,chance) # Coal ore
		repopulateSphere(level,box,options,R,(21,0),2,1,1,1,23,chance) # Lapis ore
		repopulateSphere(level,box,options,R,(56,0),1,1,1,1,12,chance) # Diamond ore
		repopulateSphere(level,box,options,R,(73,0),2,1,1,1,12,chance) # Redstone ore
		repopulateSphere(level,box,options,R,(129,0),0,0,1,1,29,chance) # Emerald ore
	elif OPER == "Ores: Ore Spheres - No Layers":
		repopulateSphere(level,box,options,R,(14,0),2,1,1,box.miny,box.maxy,chance) # Gold ore
		repopulateSphere(level,box,options,R,(15,0),2,1,1,box.miny,box.maxy,chance) # Iron ore
		repopulateSphere(level,box,options,R,(16,0),3,1,1,box.miny,box.maxy,chance) # Coal ore
		repopulateSphere(level,box,options,R,(21,0),2,1,1,box.miny,box.maxy,chance) # Lapis ore
		repopulateSphere(level,box,options,R,(56,0),1,1,1,box.miny,box.maxy,chance) # Diamond ore
		repopulateSphere(level,box,options,R,(73,0),2,1,1,box.miny,box.maxy,chance) # Redstone ore
		repopulateSphere(level,box,options,R,(129,0),0,0,1,box.miny,box.maxy,chance) # Emerald ore
	elif OPER == "Ores: Ore layers - Stained Clay":
		repopulateLayers(level,box,options,R,(159,0),R.randint(16,width),chance)
	elif OPER == "Interference":
		interference(level,box,options,R,material,R.randint(16,width),chance,False)
	elif OPER == "Interference: Variable Frequency":
		interferenceVariableFreq(level,box,options,R,material,R.randint(16,width),chance,False)
	elif OPER == "Ores: Twisting Seam":
		repopulateTwisting(level,box,options,R,(R.randint(16,width-16),R.randint(16,height-16),R.randint(16,depth-16)),40,-1,-1,material,chance)
	elif OPER == "Columns":
		columns(level,box,options,R,material,R.randint(2,int(sqrt(width**2+depth**2))))
	elif OPER == "Volcano":
		volcano(level,box,options,R,material)
	elif OPER == "Maze: Stone":
		mazeStone(level,box,options,R,material)
	elif OPER == "Tunnel: Smooth Path":
		smoothTunnel(level, box, options, R, chance)
	elif OPER == "Perlin 3D":
		perlinoise3D(level, box, options, R)
	elif OPER == "Perlin 2D":
		perlinoise2D(level, box, options, R)
	elif OPER == "Limestone Cave":
		limestoneCave(level, box, options, R)
	elif OPER == "Web":
		web(level, box, options, R)
	elif OPER == "Denseland":
		denseland(level, box, options, R)
#	print "Invokes = "+str(COUNTTUNNELINVOKES)
		
	FuncEnd(level,box,options,method)

def web(level,box,options,R):
	method = "web"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method)

	# ISSUES: Draping / curving sections to be done. Real webs aren't straight lines. Use Bresenham?
	
	# Init
	(blockID, blockData) = getBlockFromOptions(options, "Material:")
	(stickID, stickData) = getBlockFromOptions(options, "Block Type to Replace:")
	ERROR = options["Chance:"]

	SPIDERSIZE = R.randint(8,16)
	# Anchor the web to two sticks starting at these points
	(stick1x1, stick1y1, stick1z1) = (R.randint(0,centreWidth), 0, R.randint(0,centreDepth) )
	(stick2x1, stick2y1, stick2z1) = (R.randint(0,centreWidth), 0, R.randint(0,centreDepth) )
	
	# Bridge
		# Within the select box, choose two anchor points for the web
		# Ensure the points are widely spaced apart
	(x,y,z) = (ap1x, ap1y, ap1z) = (0, R.randint(centreHeight,height-1), R.randint(0,centreDepth) )
	(x1,y1,z1) = (ap2x, ap2y, ap2z) = (width-1, R.randint(centreHeight,height-1), R.randint(centreDepth,depth-1) )
	drawLine(level, (stickID, stickData), (stick1x1, stick1y1, stick1z1), (ap1x, ap1y, ap1z) )
	drawLine(level, (stickID, stickData), (stick2x1, stick2y1, stick2z1), (ap2x, ap2y, ap2z) )
	drawLine(level, (blockID, blockData), (x,y,z), (x1,y1,z1) )

	# Main spar
		# Halfway across bridge
	(x,y,z) = (ms1x,ms1y,ms1z) = ((ap1x+ap2x)/2,(ap1y+ap2y)/2,(ap1z+ap2z)/2)
		# distance along stick1
	(stick1slopex,stick1slopey,stick1slopez) = (ap1x-stick1x1,ap1y-stick1y1,ap1z-stick1z1)
	spar1maxdist = sqrt(stick1slopex**2+stick1slopey**2+stick1slopez**2)
	spar1dist = 0.0 #R.random()/2
	(x1,y1,z1) = (ap3x, ap3y, ap3z) = (stick1x1+spar1dist*stick1slopex,stick1y1+spar1dist*stick1slopey,stick1z1+spar1dist*stick1slopez)
	drawLine(level, (blockID, blockData), (x,y,z), (x1,y1,z1) )
	
	# Centre blob
	(cowx, cowy, cowz) = ((ms1x+ap3x)/2,(ms1y+ap3y)/2,(ms1z+ap3z)/2) # Centre of web	

	Q = []
	# Spars & 	# Web spiral
	for P in drawWebSpars(level,box,(blockID, blockData),SPIDERSIZE,(stick1x1, stick1y1, stick1z1),(ap1x, ap1y, ap1z),(cowx, cowy, cowz)):
		Q.append(P)
	for P in drawWebSpars(level,box,(blockID, blockData),SPIDERSIZE,(ap1x, ap1y, ap1z),(ap2x, ap2y, ap2z),(cowx, cowy, cowz)):
		Q.append(P)
	for P in drawWebSpars(level,box,(blockID, blockData),SPIDERSIZE,(ap2x, ap2y, ap2z),(stick2x1, stick2y1, stick2z1),(cowx, cowy, cowz)):
		Q.append(P)
	drawLine(level, (blockID, blockData), (stick1x1, stick1y1, stick1z1), (stick2x1, stick2y1, stick2z1) )
	for P in drawWebSpars(level,box,(blockID, blockData),SPIDERSIZE,(stick2x1, stick2y1, stick2z1),(stick1x1, stick1y1, stick1z1),(cowx, cowy, cowz)):
		Q.append(P)
	# Now I have an array of points in space. Join them up!
	R = []
	for T in Q:
		S = []
		#print T
		for i in xrange(0,len(T)):
			if len(R) > i:
				drawLine(level, (stickID, stickData), T[i], R[(i)%len(R)] )
			#print T[i]
			S.append(T[i])
		R = S
#		for (x,y,z) in S:
#			R.append((x,y,z))
#		print R
	

	# Centre eat and mesh
	
	FuncEnd(level,box,options,method)	

def drawWebSpars(level,box,material,gap,(p1x,p1y,p1z),(p2x,p2y,p2z),(cowx,cowy,cowz)):
	(stick1slopex,stick1slopey,stick1slopez) = (p2x-p1x,p2y-p1y,p2z-p1z)
	spar1maxdist = sqrt(stick1slopex**2+stick1slopey**2+stick1slopez**2)
	spar1dist = 0.0 #R.random()/2
	# Remaining spars
	unit = 1.0/spar1maxdist
	keepGoing = True
	count = 0
	Q = []
	while keepGoing == True and count < 10000:
		count = count +1 # Limit looping so we don't run forever
		if count%10 == 0:
			print count
		# Based on the spider's reach, build the spars from the centre outward.
		# Traverse Stick 1 at nominated intervals, draw to the middle: (ms1x,ms1y,ms1z)
		spar1dist = spar1dist+gap*unit
		P = []
		if spar1dist <= 1:
			(x1,y1,z1) = (p1x+spar1dist*stick1slopex,p1y+spar1dist*stick1slopey,p1z+spar1dist*stick1slopez)
			drawLine(level, material, (x1,y1,z1), (cowx,cowy,cowz) )
			(dx,dy,dz) = (x1-cowx,y1-cowy,z1-cowz)
			d = sqrt(dx**2+dy**2+dz**2)
			du = 1.0/d
			ds = 0.0
			kg = True
			while kg == True:
				ds = ds + gap*du # Try gap/3 for fun!
				if ds <= 1.0:
					P.append( (cowx+dx*ds,cowy+dy*ds,cowz+dz*ds) ) # add each point along the spar at the gap
				else:
					kg = False				
			Q.append(P) # Add in the set of points for this strut
		else:
			keepGoing = False
	
	return Q
	
def getBlockFromOptions(options, block):
	return (options[block].ID, options[block].blockData)

def drawLine(scratchpad, (blockID, blockData), (x,y,z), (x1,y1,z1) ):
	drawLineConstrained(scratchpad, (blockID, blockData), (x,y,z), (x1,y1,z1), 0 )

def drawLine1(scratchpad, (blockID, blockData), (x,y,z), (x1,y1,z1) ):
	for px, py, pz in bresenham.bresenham((x,y,z),(x1,y1,z1)):
		setBlock(scratchpad,(blockID, blockData),px,py,pz)
	setBlock(scratchpad,(blockID, blockData),x1,y1,z1)

def drawLineConstrained(scratchpad, (blockID, blockData), (x,y,z), (x1,y1,z1), maxLength ):
	dx = x1 - x
	dy = y1 - y
	dz = z1 - z

	distHoriz = dx*dx + dz*dz
	distance = sqrt(dy*dy + distHoriz)

	if distance < maxLength or maxLength < 1:
		phi = atan2(dy, sqrt(distHoriz))
		theta = atan2(dz, dx)

		iter = 0
		while iter <= distance:
			scratchpad.setBlockAt((int)(x+iter*cos(theta)*cos(phi)), (int)(y+iter*sin(phi)), (int)(z+iter*sin(theta)*cos(phi)), blockID)
			scratchpad.setBlockDataAt((int)(x+iter*cos(theta)*cos(phi)), (int)(y+iter*sin(phi)), (int)(z+iter*sin(theta)*cos(phi)), blockData)
			iter = iter+0.5 # slightly oversample because I lack faith.
	
def denseland(level,box,options,R):
	# Ok this isn't really Perlin. It's my quick kludge
	method = "denseland"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method)

	maxDIST = float((centreWidth+centreDepth+centreHeight)/3)**2
	
	noiseAbove3d = makeSomeNoiseRange(width, height, depth, R, -1,1)
	spaceAbove3d = zeros((width,height,depth)) # density field in 3D space
	perlman(noiseAbove3d,spaceAbove3d,[2,4],2.0)
	
	noise3d = makeSomeNoiseRange(width, height, depth, R, -1,1)
	space3d = zeros((width,height,depth)) # density field in 3D space
	perlman(noise3d,space3d,[2,4],2.0)

	surfaceNoise2d = makeSomeNoiseRange(width, 1, depth, R, -1,1)
	surface2d = zeros((width,1,depth)) # density field in 2D space
	perlman(surfaceNoise2d,surface2d,[16,32],sqrt(2.0))

	
	# apply any filters
	for x in xrange(0,width):
		if x%10 == 0:
			print x
		for z in xrange(0,depth):
			for y in xrange(0,height):
				dist = (x-centreWidth)**2+(z-centreDepth)**2+(y-centreHeight)**2
				
				if dist <= maxDIST:
					space3d[x][y][z] = space3d[x][y][z]+(-1.0+1*dist/maxDIST)

				
				
				
	# merge
	newSpace = zeros((width,height,depth))
	for x in xrange(0,width):
		if x%10 == 0:
			print x
		for z in xrange(0,depth):
			ySurf = int(centreHeight+float(surface2d[x][0][z]*centreHeight))%height
			#print ySurf
			
			for y in xrange(0,height):
				if y < ySurf:
					newSpace[x][y][z] = space3d[x][y][z]
				elif y == ySurf and space3d[x][y][z] > -0.8:
					newSpace[x][y][z] = 10
				# elif y-ySurf > 32:
#					if spaceAbove3d[x][y-ySurf][z] > 0.8:
#						newSpace[x][y][z] = spaceAbove3d[x][y-ySurf][z] #-1.3*sin(pi/2*(y-ySurf)/16) # crash if ySurf too high
#					else:
#						newSpace[x][y][z] = -1000
				else:
					newSpace[x][y][z] = -1000
#				if y == height -1:	newSpace[x][y][z] = 0.9 # Cloud layer goes here
	# draw	
	drawSpace(level,box,options,newSpace)
	FuncEnd(level,box,options,method)	

def drawSpace(level,box,options,space):				
	method = "drawSpace"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method)
	(mx,my,mz) = space.shape
				
	for x in xrange(0,width):
		if x%10 == 0:
			print x
		for z in xrange(0,depth):
			for y in xrange(0,height):				
				n = space[x%mx][y%my][z%mz]
				if n >= 10:
					setBlock(level,(2,0),box.minx+x,box.miny+y,box.minz+z) # Grass
				elif n > 0.8:
					setBlock(level,(155,0),box.minx+x,box.miny+y,box.minz+z)
				elif n > 0.6:
					setBlock(level,(1,6),box.minx+x,box.miny+y,box.minz+z)
				elif n > 0.5:
					setBlock(level,(1,1),box.minx+x,box.miny+y,box.minz+z)
				elif n > 0.0:
					setBlock(level,(1,2),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.2:
					setBlock(level,(1,3),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.3:
					setBlock(level,(1,4),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.4:
					setBlock(level,(1,5),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.8:
					setBlock(level,(1,0),box.minx+x,box.miny+y,box.minz+z)
	FuncEnd(level,box,options,method)	
	
def perlman(noise,space,octaves,persistence):
	# Ok this isn't really Perlin. It's my quick kludge
	method = "perlmanoise"
	print method+' started'
	(mx,my,mz) = space.shape
	(width,height,depth) = (mx,my,mz)
	(centreWidth,centreDepth,centreHeight) = (int(mx/2),int(my/2),int(mz/2))
	if persistence <=0:
		persistence = 2.0
	if len(octaves) < 1:
		octaves = [2,4]
	maxDIST = float((centreWidth+centreDepth+centreHeight)/3)**2
	
	for x in xrange(0,mx):
		if x%10 == 0:
			print x
		for z in xrange(0,mz):
			for y in xrange(0,my):
				dist = (x-centreWidth)**2+(z-centreDepth)**2+(y-centreHeight)**2
				n=0.0
				for octave in octaves:
					stepSize = persistence*octave
					px_int = int(x/stepSize) 	# Noise sample location
					pz_int = int(z/stepSize)
					py_int = int(y/stepSize)
					pdx_float = x-(px_int*stepSize)	# position within the cell
					pdz_float = z-(pz_int*stepSize)
					pdy_float = y-(py_int*stepSize)
					# I need to find what the noise value is within this cell.
					rdx = pdx_float/stepSize #0.5+cos(pi+pi*dx/v_int)/2
					rdz = pdz_float/stepSize #0.5+cos(pi+pi*dz/v_int)/2
					rdy = pdy_float/stepSize #0.5+cos(pi+pi*dz/v_int)/2
					
					p000 = noise[px_int][py_int][pz_int]
					p001 = noise[px_int][py_int][(pz_int+1)%depth]
					p010 = noise[px_int][(py_int+1)%height][pz_int]
					p011 = noise[px_int][(py_int+1)%height][(pz_int+1)%depth]					
					p100 = noise[(px_int+1)%width][py_int][pz_int]					
					p101 = noise[(px_int+1)%width][py_int][(pz_int+1)%depth]					
					p110 = noise[(px_int+1)%width][(py_int+1)%height][pz_int]					
					p111 = noise[(px_int+1)%width][(py_int+1)%height][(pz_int+1)%depth]					
			
					dvx1 = cosineInterpolate(p000,p100,rdx) #
					dvx2 = cosineInterpolate(p010,p110,rdx)
					dvx3 = cosineInterpolate(p001,p101,rdx) #
					dvx4 = cosineInterpolate(p011,p111,rdx)

					dvz1 = cosineInterpolate(dvx1,dvx3,rdz) #
					dvz2 = cosineInterpolate(dvx2,dvx4,rdz) #
					
					if dist <= maxDIST:
						n = n+cosineInterpolate(dvz1,dvz2,rdy)	
					else:
						n = n+cosineInterpolate(dvz1,dvz2,rdy)
					
				# Save it
				n = n/len(octaves)
				space[x][y][z] = n

	print method+' ended'
	
def limestoneCave(level,box,options,R):
	# Ok this isn't really Perlin. It's my quick kludge
	method = "limestoneCave"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method)	
	SQRT2 = sqrt(2)
	# Raise the roof (Make some noise)
	noise = makeSomeNoiseRange(width, height, depth, R, -1,1)
	persistence = 2.0
	L = [2,4]
	maxDIST = ((centreWidth+centreDepth+centreHeight)/3)**2

	for x in xrange(0,width):
		if x%10 == 0:
			print x
		for z in xrange(0,depth):
			for y in xrange(0,height):
				dist = (x-centreWidth)**2+(z-centreDepth)**2+(y-centreHeight)**2
				n=0.0
				for octave in L:
					stepSize = persistence*octave
					px_int = int(x/stepSize) 	# Noise sample location
					pz_int = int(z/stepSize)
					py_int = int(y/stepSize)
					pdx_float = x-(px_int*stepSize)	# position within the cell
					pdz_float = z-(pz_int*stepSize)
					pdy_float = y-(py_int*stepSize)
					# I need to find what the noise value is within this cell.
					rdx = pdx_float/stepSize #0.5+cos(pi+pi*dx/v_int)/2
					rdz = pdz_float/stepSize #0.5+cos(pi+pi*dz/v_int)/2
					rdy = pdy_float/stepSize #0.5+cos(pi+pi*dz/v_int)/2
					
					p000 = noise[px_int][py_int][pz_int]
					p001 = noise[px_int][py_int][(pz_int+1)%depth]
					p010 = noise[px_int][(py_int+1)%height][pz_int]
					p011 = noise[px_int][(py_int+1)%height][(pz_int+1)%depth]					
					p100 = noise[(px_int+1)%width][py_int][pz_int]					
					p101 = noise[(px_int+1)%width][py_int][(pz_int+1)%depth]					
					p110 = noise[(px_int+1)%width][(py_int+1)%height][pz_int]					
					p111 = noise[(px_int+1)%width][(py_int+1)%height][(pz_int+1)%depth]					
			
					dvx1 = linearInterpolate(p000,p100,rdx) #
					dvx2 = linearInterpolate(p010,p110,rdx)
					dvx3 = linearInterpolate(p001,p101,rdx) #
					dvx4 = linearInterpolate(p011,p111,rdx)

					dvz1 = linearInterpolate(dvx1,dvx3,rdz) #
					dvz2 = linearInterpolate(dvx2,dvx4,rdz) #
					
					if dist <= maxDIST:
						n = n+cosineInterpolate(dvz1,dvz2,rdy)+(-2.0+2*dist/maxDIST)	
					else:
						n = n+cosineInterpolate(dvz1,dvz2,rdy)
					
				# Draw it
				n = n/len(L)
				#if y < centreHeight:
				if n > 0.8:
					setBlock(level,(155,0),box.minx+x,box.miny+y,box.minz+z)
				elif n > 0.6:
					setBlock(level,(1,0),box.minx+x,box.miny+y,box.minz+z)
				elif n > 0.5:
					setBlock(level,(1,1),box.minx+x,box.miny+y,box.minz+z)
				elif n > 0.0:
					setBlock(level,(1,2),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.2:
					setBlock(level,(1,3),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.3:
					setBlock(level,(1,4),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.4:
					setBlock(level,(1,5),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.41:
					setBlock(level,(89,0),box.minx+x,box.miny+y,box.minz+z)				
#				elif n > -0.5:
#					setBlock(level,(20,0),box.minx+x,box.miny+y,box.minz+z)
				else:
					setBlock(level,(0,0),box.minx+x,box.miny+y,box.minz+z)
				#else:
				#	if n < -0.8:
				#		setBlock(level,(1,3),box.minx+x,box.miny+y,box.minz+z)
					
	FuncEnd(level,box,options,method)

def linearInterpolate(a, b, x):
	return a+(b-a)*x

def cosineInterpolate(a, b, x): # http://www.minecraftforum.net/forums/off-topic/computer-science-and-technology/482027-generating-perlin-noise?page=40000
	ft = pi * x
	f = ((1.0 - cos(ft)) * 0.5)
	ret = a * (1.0 - f) + b * f
	return ret

def perlinoise2D(level,box,options,R):
	# Ok this isn't really Perlin. It's my quick kludge
	method = "perlmanoise2D"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method)	
	SQRT2 = sqrt(2)
	# Raise the roof (Make some noise)
	noise = makeSomeNoiseRange(width, 1, depth, R, -1,1)
	#noise = smoothNoise2D(noise,3)
	#noise = smoothNoise2D(noise,1)
	persistence = 1.41
	L = [1,8,16,32,64]
	for x in xrange(0,width):
		if x%10 == 0:
			print x
		for z in xrange(0,depth):
			n=0.0
			stepSizeMin = 1.0
			for octave in L:
				v = persistence*octave
				v_int = int(v)
				stepSize = persistence*octave
				
				px_int = int(x/stepSize) 	# Noise sample location
				pz_int = int(z/stepSize)
				pdx_float = x-(px_int*stepSize)	# position within the cell
				pdz_float = z-(pz_int*stepSize)
				# I need to find what the noise value is within this cell.
				rdx = pdx_float/stepSize #0.5+cos(pi+pi*dx/v_int)/2
				rdz = pdz_float/stepSize #0.5+cos(pi+pi*dz/v_int)/2
				
				p00 = noise[px_int][0][pz_int]
				p10 = noise[(px_int+1)%width][0][pz_int]
				p01 = noise[px_int][0][(pz_int+1)%depth]
				p11 = noise[(px_int+1)%width][0][(pz_int+1)%depth]

				dvx1 = cosineInterpolate(p00,p10,rdx)
				dvx2 = cosineInterpolate(p01,p11,rdx)

				n = n+cosineInterpolate(dvx1,dvx2,rdz)*float(octave)

			maxV = int(centreHeight+n)
			for i in xrange(int(maxV)-16,int(maxV)):
				material = (1,0)
				if i >= maxV-1:
					material =  (155,0) # blocks[int(value*len(blocks)%len(blocks))]
		
				setBlock(level,material,box.minx+x,box.miny+i,box.minz+z)
			if int(maxV) < centreHeight:
				for i in xrange(int(maxV),int(centreHeight)):
					setBlock(level,(9,0),box.minx+x,box.miny+i,box.minz+z)
		
	FuncEnd(level,box,options,method)

def smooshNoise3D(noise,ythresh,val):
	(mx,my,mz) = noise.shape
	smooshNoise = zeros((noise.shape))
	
	for x in xrange(0,mx):
		for z in xrange(0,mz):
			for y in xrange(0,my):
				smooshNoise[x][y][z] = noise[x][y][z]
				if y > ythresh:
					smooshNoise[x][y][z] = smooshNoise[x][y][z]/val
	return smooshNoise	
	
def perlinoise3D(level,box,options,R):
	# Ok this isn't really Perlin. It's my quick kludge
	method = "perlmanoise3D"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method)	
	SQRT2 = sqrt(2)
	# Raise the roof (Make some noise)
	noise = makeSomeNoiseRange(width, height, depth, R, -1,1)
	noise = smooshNoise3D(noise,int(height/2),10)
	#print noise
	#noise = smoothNoise2D(noise,3)
	#noise = smoothNoise2D(noise,1)
	persistence = 2.0
	L = [2,4]
	for x in xrange(0,width):
		if x%10 == 0:
			print x
		for z in xrange(0,depth):
			for y in xrange(0,height):
				n=0.0
				stepSizeMin = 1.0
				for octave in L:
					#v = persistence**octave
					#v_int = int(v)
					stepSize = persistence*octave
					px_int = int(x/stepSize) 	# Noise sample location
					pz_int = int(z/stepSize)
					py_int = int(y/stepSize)
					pdx_float = x-(px_int*stepSize)	# position within the cell
					pdz_float = z-(pz_int*stepSize)
					pdy_float = y-(py_int*stepSize)
					# I need to find what the noise value is within this cell.
					rdx = pdx_float/stepSize #0.5+cos(pi+pi*dx/v_int)/2
					rdz = pdz_float/stepSize #0.5+cos(pi+pi*dz/v_int)/2
					rdy = pdy_float/stepSize #0.5+cos(pi+pi*dz/v_int)/2
					
					p000 = noise[px_int][py_int][pz_int]
					p001 = noise[px_int][py_int][(pz_int+1)%depth]
					p010 = noise[px_int][(py_int+1)%height][pz_int]
					p011 = noise[px_int][(py_int+1)%height][(pz_int+1)%depth]					
					p100 = noise[(px_int+1)%width][py_int][pz_int]					
					p101 = noise[(px_int+1)%width][py_int][(pz_int+1)%depth]					
					p110 = noise[(px_int+1)%width][(py_int+1)%height][pz_int]					
					p111 = noise[(px_int+1)%width][(py_int+1)%height][(pz_int+1)%depth]					
			
					# print str(x)+' '+str(y)+' '+str(z)+' '+str(p000)+' '+str(p001)+' '+str(p010)+' '+str(p011)+' '+str(p100)+' '+str(p101)+' '+str(p110)+' '+str(p111)
			
					dvx1 = cosineInterpolate(p000,p100,rdx) #
					dvx2 = cosineInterpolate(p010,p110,rdx)
					dvx3 = cosineInterpolate(p001,p101,rdx) #
					dvx4 = cosineInterpolate(p011,p111,rdx)

					#dvy1 = cosineInterpolate(p000,p010,rdy)
					#dvy2 = cosineInterpolate(p100,p110,rdy)
					#dvy3 = cosineInterpolate(p001,p011,rdy)
					#dvy4 = cosineInterpolate(p101,p111,rdy)

					dvz1 = cosineInterpolate(dvx1,dvx3,rdz) #
					dvz2 = cosineInterpolate(dvx2,dvx4,rdz) #
					#dvz3 = cosineInterpolate(p010,p011,rdz)
					#dvz4 = cosineInterpolate(p110,p111,rdz)
					
					dist = sqrt((centreWidth-x)**2+(centreDepth-z)**2+(centreHeight-y)**2)
					maxDIST = sqrt((centreWidth)**2+(centreDepth)**2+(centreHeight)**2)
					n = n+cosineInterpolate(dvz1,dvz2,rdy)-(sin(pi*y/height))#+0.5*cos(pi*dist/maxDIST)
					#(2.0*sin(y/height*pi/2))
					#print str(dvx1)+' '+str(dvx2)
					#dvz1 = noise[px_int][0][(pz_int+1)%depth]-noise[px_int][0][pz_int]
					#dvz2 = noise[(px_int+1)%width][0][(pz_int+1)%depth]-noise[(px_int+1)%width][0][pz_int]
					# Now we have the average 'delta' vectors in the x and z direction for this noise cell
					
					# Where are we in the cell, and what is the noise value here?
					#nx = (rdx*dvx1+(1.0-rdx)*dvx2)/2+noise[px_int][0][pz_int] # Funny twisty surface segment - linear interpolation with a twisty bit
					#nz = (rdz*dvz1+(1.0-rdz)*dvz2)/2+noise[px_int][0][pz_int]
					#n = (nx+nz)/2
					#print noise[(px_int+1)%width][0][(pz_int+1)%depth]
					#print noise[(px_int+1)%width][0][pz_int]
					#print n
					# Draw it
				n = n/len(L)
				#if y < centreHeight:
				if n > 0.8:
					setBlock(level,(155,0),box.minx+x,box.miny+y,box.minz+z)
				elif n > 0.6:
					setBlock(level,(1,0),box.minx+x,box.miny+y,box.minz+z)
				elif n > 0.5:
					setBlock(level,(1,1),box.minx+x,box.miny+y,box.minz+z)
				elif n > 0.0:
					setBlock(level,(1,2),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.2:
					setBlock(level,(1,3),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.3:
					setBlock(level,(1,4),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.4:
					setBlock(level,(1,5),box.minx+x,box.miny+y,box.minz+z)
				elif n > -0.8:
					setBlock(level,(1,6),box.minx+x,box.miny+y,box.minz+z)
				#else:
				#	if n < -0.8:
				#		setBlock(level,(1,3),box.minx+x,box.miny+y,box.minz+z)
					
	FuncEnd(level,box,options,method)
	
def makeSomeNoiseRange(width, height, depth, R, a, b): # Other function names considered were: ShoutOutLoud and BoomBoomShakeTheRoom
	method = "makeSomeNoise"
	# create a noise 
	noise = zeros((width,height,depth))
	for x in range(width):
		for z in range(depth):
			for y in range(height):
				noise[x][y][z] = a+(b-a)*R.random()
	#			print noise[x][y][z]
	return noise
	
def makeSomeNoise(width, height, depth, R): # Other function names considered were: ShoutOutLoud and BoomBoomShakeTheRoom
	method = "makeSomeNoise"
	# create a noise 
	noise = zeros((width,height,depth))
	for x in range(width):
		for z in range(depth):
			for y in range(height):
				noise[x][y][z] = R.random()
	return noise

def smoothNoise2D(noise,dist):
	(mx,my,mz) = noise.shape
	
	smoothNoise = zeros((noise.shape))
	
	for x in xrange(0,mx):
		for z in xrange(0,mz):
			for y in xrange(0,my):
				valhere = 0.0
				numPts = 0
				for r in xrange(-dist,dist):
					for s in xrange(0,1):
						for t in xrange(-dist,dist):
							valhere = valhere+noise[(x+r)%mx][(y+s)%my][(z+t)%mz]
							numPts += 1
				valhere = valhere / numPts
				smoothNoise[x][y][z] = valhere
	return smoothNoise
	
def smoothTunnel(level, box, options, R, chance):
	method = "smoothTunnel"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method)	

	F1 = 0.25
	F2 = 0.75
	DRAWSOURCE = False
	MAXSIZE = 2
	sizeX = MAXSIZE
	sizeY = MAXSIZE
	sizeZ = MAXSIZE

	
	P = []
	count = 0
	(xo,yo,zo) = (-1,-1,-1)
	(x0,y0,z0) = (-1,-1,-1)
	for i in xrange(0,R.randint(7,170)): # 

		(x1,y1,z1) = (R.randint(0,width),R.randint(0,height),R.randint(0,depth))
		if count == 0:
			(xo,yo,zo) = (x1,y1,z1)
		count = count+1
		P.append( (x1,y1,z1) )
		(x0,y0,z0) = (x1,y1,z1)
	P.append( (xo,yo,zo) ) # loop around to the first point

		
	for i in xrange(0,5):
		P = chaikinSmoothAlgorithm(P)
		
	(x0,y0,z0) = (-1,-1,-1)
	count = 0
	for (px,py,pz) in P:
		if count > 0: # We have a previous point
			for y in xrange(-sizeY,sizeY+1):
				for x in xrange(-sizeX,sizeX+1):
					for z in xrange(-sizeZ,sizeZ+1):
						if not ((abs(x) == sizeX and abs(y) == sizeY) or (abs(y) == sizeY and abs(z) == sizeZ) or (abs(x) == sizeX and abs(z) == sizeZ)): # taper walls
							drawLine(level,AIR,(box.minx+x0+x,box.miny+y0+y,box.minz+z0+z),(box.minx+px+x,box.miny+py+y,box.minz+pz+z))
		else:
			count = count+1
		(x0,y0,z0) = (px,py,pz)
	
	FuncEnd(level,box,options,method)

def sphere(level,box,options,material,(px,py,pz),radius):
	method = "Sphere"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start

	for y in xrange(-radius,radius+1):
		for x in xrange(-radius,radius+1):
			for z in xrange(-radius,radius+1):
				if y**2+z**2+x**2 <= radius**2:
					setBlock(level,material,box.minx+px+x,box.miny+py+y,box.minz+pz+z)

	FuncEnd(level,box,options,method)		

def tunnelSection(level,box,options,R,(px,py,pz),radius,MAXLEN,theta,phi,chance,wiggle,edge):
	method = "Tunnel:Section"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start

	MAXSIZE = radius
	sizeX = MAXSIZE
	sizeY = MAXSIZE
	sizeZ = MAXSIZE
	velocity = 3	
	steps = 0
	MAXSTEPS = MAXLEN
	if MAXLEN < 1:
		MAXSTEPS = R.randint(10,50)

	keepGoing = True
	while keepGoing == True:
		if R.randint(1,100) < chance:
			# Draw tunnel segment
			for y in xrange(-sizeY,sizeY+1):
				for x in xrange(-sizeX,sizeX+1):
					for z in xrange(-sizeZ,sizeZ+1):
						if not ((abs(x) == sizeX and abs(y) == sizeY) or (abs(y) == sizeY and abs(z) == sizeZ) or (abs(x) == sizeX and abs(z) == sizeZ)): # taper walls
							setBlock(level,AIR,box.minx+px+x,box.miny+py+y,box.minz+pz+z)	

		if wiggle == True:
			# Change direction and speed slightly
			theta = theta+R.randint(-3,3)*a*20 # Sharp turns wrap it up tightly
			if R.randint(0,100) > 50:
				phi = phi+R.randint(-1,1)*a*2							

		velocity = velocity+R.randint(-1,2)
		if velocity < 2:
			velocity = 2
		elif velocity > 3:
			velocity = 3

		#  Move to new segment, repeat
		px = px+velocity*cos(theta)*cos(phi)
		py = py+velocity*sin(phi)
		pz = pz+velocity*sin(theta)*cos(phi)

		if edge == True:
			# Bounds checking
			if px < 16 or pz < 16 or px >= width-16 or pz >= depth-16:
				theta = theta + pi
			if py < 16 or py > height-1:
				phi = phi + pi
		
		steps = steps+1
		if steps > MAXSTEPS:
			keepGoing = False

	return (px,py,pz)
		
def volcano(level,box,options,R,material): # See https://en.wikipedia.org/wiki/Volcano for a diagram
	method = "Volcano"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start

	# Make a cinder cone
	(px,py,pz) = (centreWidth,0,centreDepth)
	maxheight = height-1
	y = maxheight # R.randint(int(maxheight/3),maxheight)
	if y > centreWidth:
		y = int(centreWidth-1)
	if y > centreDepth:
		y = int(centreDepth-1)
	maxheight = y
	
	base = 1
	keepGoing = True
	skip = R.randint(int(maxheight/16),int(maxheight/4)) # This is the number of layers not to draw at the top of the cone.
	while keepGoing == True and y > 0:
		if y%10 == 0:
			print y
		y = y - 1
		skip = skip -1
		if skip < 1:
			aircount = 0
			r = maxheight-y
			# Draw a filled disc at this layer of radius r
			for x in xrange(-r,r):
				for z in xrange(-r,r):
					r1 = x**2+z**2
					if sqrt(r1) <= r:
						b = getBlock(level, box.minx+px+x, box.miny+py+y, box.minz+pz+z)
						#print str(box.minx)+','+str(box.miny)+','+str(box.minz)
						#print 'Block:'+str(b)+' '+str(r)+' '+str(box.minx+px+x)+','+str(box.miny+py+y)+','+str(box.minz+pz+z)
						if b == AIR: # air there
							aircount = aircount+1
							setBlock(level, material, box.minx+px+x, box.miny+py+y, box.minz+pz+z)
							#print aircount
		
			# If there was nothing to draw this time, we've reached the base (surface) and can stop
			if aircount  == 0:
				keepGoing = False
				base = y

	(px,py,pz) = (centreWidth,R.randint(int(base/3),base),centreDepth)
			
	# Make the conduit and throat
	theta = R.randint(0,359)
	phi = pi/2 # Straight up, dawg.
	tunnelSection(level,box,options,R,(px,py,pz),R.randint(2,4),height,theta,phi,100,False,False)
	
	# Throat
	offset = int((maxheight - py)/2)
	for i in xrange(3,11):
		theta = R.randint(0,359)
		phi = pi/2-randint(1,10)*a # Straight up, dawg.
		tunnelSection(level,box,options,R,(px,py+offset,pz),R.randint(2,4),height,theta,phi,100,False,False)
		
	# Make Dikes
	print "Dikes"
	for i in xrange(1,6):
		theta = R.randint(0,359)
		phi = R.random()*pi/2+pi/2
		tunnelSection(level,box,options,R,(px,py,pz),2,centreHeight,theta,phi,100,True,False)

	# Make a magma chamber
	print "Magma chamber"
	for i in xrange(3,R.randint(6,15)):
		d = 1+i
		(dx,dy,dz) = (R.randint(-d,d),R.randint(-d,d),R.randint(-d,d))
		r = R.randint(2*i,3*i)
		sphere(level,box,options,LAVA,(px+dx,py+dy,pz+dz),r)

	
	
	FuncEnd(level,box,options,method)		

def columns(level,box,options,R,material,count):
	method = "Columns"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start

	for iter in xrange(1,count):
		# print iter
		if width >0 and depth > 0:
			column(level,box,options,R,material,(R.randint(0,width),R.randint(0,depth)),R.randint(1,2))

	FuncEnd(level,box,options,method)		
		
def column(level,box,options,R,material,(px,pz),radius):
	method = "Column"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start
	
	dxm = R.randint(1,radius)
	dxp = R.randint(1,radius)
	dzm = R.randint(1,radius)
	dzp = R.randint(1,radius)
	
	for x in xrange(-dxm,dxp):
		for z in xrange(-dzm,dzp):
			r = x**2+z**2

			if sqrt(r) <= radius:
				mh = abs(int(height/2)-r)
				
				if getBlock(level,box.minx+px+x,box.miny-1,box.minz+pz+z) != AIR:
					for y in xrange(0,randint(0,mh)):
						setBlock(level,material,box.minx+px+x,box.miny+y,box.minz+pz+z)	
				if getBlock(level,box.minx+px+x,box.maxy,box.minz+pz+z) != AIR:
					for y in xrange(0,randint(0,mh)):
						setBlock(level,material,box.minx+px+x,box.maxy-y,box.minz+pz+z)	
		
	FuncEnd(level,box,options,method)		
	
def repopulateLayers(level,box,options,R,(blockID,blockData),wavelen,chance):
	method = "Repopulate:Layers"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start
	
	P = []
	for i in xrange(0, R.randint(2,5)):
		p = (R.randint(0,width),R.randint(0,height),R.randint(0,depth))
		P.append(p)
	
	replaceMaterialID = options["Block Type to Replace:"].ID
	
	for x in xrange(0,width):
		if x%10 == 0:
			print x
		for z in xrange(0,depth):
			for y in xrange(0,height):
				block = level.blockAt(box.minx+x,box.miny+y,box.minz+z)
				if block == replaceMaterialID and R.randint(0,99) < chance: # Stone
					amp = float(0)
					for (px,py,pz) in P:
						(dx,dy,dz) = (x-px,y-py,z-pz)
						dist = sqrt(dx**2+dy**2+dz**2)
						amp = amp + cos(dist*a/(2*pi)*wavelen)
					amp = (abs(15*amp))/len(P)
					setBlock(level,(blockID,int(amp)),box.minx+x,box.miny+y,box.minz+z)	

	FuncEnd(level,box,options,method)		

def interference(level,box,options,R,(blockID,blockData),wavelen,chance,REPLACE):
	method = "Interference"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start

	blocks = palette[R.randint(0,len(palette)-1)]
	DW = width+1000
	DH = height+1000
	DD = depth+1000
	P = []
	for i in xrange(0, R.randint(2,5)):
		if R.randint(0,100) > 95:
			p = (R.randint(0,width),R.randint(0,height),R.randint(0,depth))
			P.append(p)
		p = (R.randint(-DW,DW),R.randint(-DH,DH),R.randint(-DD,DD))
		P.append(p)
	
	replaceMaterialID = options["Block Type to Replace:"].ID
	
	for x in xrange(0,width):
		if x%10 == 0:
			print x
		for z in xrange(0,depth):
			for y in xrange(0,height):
				block = level.blockAt(box.minx+x,box.miny+y,box.minz+z)
				if REPLACE == False or block == replaceMaterialID and R.randint(0,99) < chance: # Stone
					amp = float(0)
					for (px,py,pz) in P:
						(dx,dy,dz) = (x-px,y-py,z-pz)
						dist = sqrt(dx**2+dy**2+dz**2)
						amp = amp + cos(dist*a/(2*pi)*wavelen)
					amp = (abs((len(colours)-1)*amp))/len(P)
					setBlock(level,blocks[int(amp)%len(blocks)],box.minx+x,box.miny+y,box.minz+z)	# (blockID,colours[int(amp)%len(colours)])

	FuncEnd(level,box,options,method)	

def interferenceVariableFreq(level,box,options,R,(blockID,blockData),wavelen,chance,REPLACE):
	method = "InterferenceVariablefreq"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start

	ferment = True
	if R.randint(0,100) > 50:
		ferment = False
	
	blocks = palette[R.randint(0,len(palette)-1)]
	DW = width+1000
	DH = height+1000
	DD = depth+1000
	A = 8

	P = []
	maxFreq = 0
	for i in xrange(0, R.randint(2,5)):
		if R.randint(0,100) > 95:
			freq = R.randint(2,A)
			maxFreq += freq
			p = (R.randint(0,width),R.randint(0,height),R.randint(0,depth),freq)
			P.append(p)
		freq = R.randint(2,A)
		maxFreq += freq
		p = (R.randint(-DW,DW),R.randint(-DH,DH),R.randint(-DD,DD),freq)
		P.append(p)
	
	replaceMaterialID = options["Block Type to Replace:"].ID
	
	for x in xrange(0,width):
		if x%10 == 0:
			print x
		for z in xrange(0,depth):
			for y in xrange(0,height):
				block = level.blockAt(box.minx+x,box.miny+y,box.minz+z)
				if REPLACE == False or block == replaceMaterialID and R.randint(0,99) < chance: # Stone
					amp = float(0)
					for (px,py,pz,f) in P:
						(dx,dy,dz) = (x-px,y-py,z-pz)
						dist = sqrt(dx**2+dy**2+dz**2)
						amp = amp + f*cos(dist*a/(2*pi)*wavelen)
					amp = (abs((len(colours)-1)*(amp/maxFreq)))/len(P) #(len(P)*maxFreq)
					col = int(amp)%len(blocks)
					if ferment == False or (ferment == True and col <= 1):
						setBlock(level,blocks[col],box.minx+x,box.miny+y,box.minz+z)	# (blockID,colours[int(amp)%len(colours)])

	FuncEnd(level,box,options,method)	
	
def repopulateSphere(level,box,options,R,material,maxradius,mincount,maxcount,minHeight,maxHeight,chance):
	method = "Repopulate:Sphere"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start
	
	# For each chunk, 
	# for (chunk, slices, point) in level.getChunkSlices(box):
		# find the centre position
	
	replaceMaterialID = options["Block Type to Replace:"].ID
	
	count = R.randint(mincount,maxcount)
	if count > 0:
		for i in xrange(0,count):
			size = R.randint(0,maxradius)
			for (cx,cz) in box.chunkPositions:
				cx = cx * CHUNKSIZE
				cz = cz * CHUNKSIZE
				min = minHeight+maxradius
				max = maxHeight-maxradius
				if max > min:
					(px,py,pz) = (R.randint(0,CHUNKSIZE-1),
									R.randint(minHeight+maxradius,maxHeight-maxradius),
									R.randint(0,CHUNKSIZE-1))
					print px,py,pz
					for y in xrange(-size,size+1):
						for x in xrange(-size,size+1):
							for z in xrange(-size,size+1):
								if y**2+z**2+x**2 <= size**2 and R.randint(0,99) < chance:
									block = level.blockAt(box.minx+px+cx+x,box.miny+py+y,box.minz+pz+cz+z)
									#if block == 1: # Stone
									if block == replaceMaterialID: # AIR - debug
										setBlock(level,material,box.minx+px+cx+x,box.miny+py+y,box.minz+pz+cz+z)	
		
	
	FuncEnd(level,box,options,method)	

def repopulateTwisting(level,box,options,R,(px,py,pz),MAXLEN,theta,phi,material,chance):
	method = "Repopulate:Twisting"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start

	PLUNGE = 32
	
	replaceMaterialID = options["Block Type to Replace:"].ID
	
	# Direction facing
	if theta == -1:
		theta = a*R.randint(0,MAXANGLES-1)
	if phi == -1:
		phi = a*R.randint(0,MAXANGLES-1)
	velocity = 3
	
	MAXSIZE = 2
	sizeX = 2
	sizeY = 2
	sizeZ = 2

	steps = 0
	MAXSTEPS = MAXLEN
	if MAXLEN < 1:
		MAXSTEPS = R.randint(10,100) # default for joker input
	keepGoing = True
	while keepGoing == True:
		if R.randint(1,100) > 1:
			# Draw tunnel segment
			for y in xrange(-sizeY,sizeY+1):
				for x in xrange(-sizeX,sizeX+1):
					for z in xrange(-sizeZ,sizeZ+1):
						#if y**2+z**2+x**2 <= 3**2:
						if not ((abs(x) == sizeX and abs(y) == sizeY) or (abs(y) == sizeY and abs(z) == sizeZ) or (abs(x) == sizeX and abs(z) == sizeZ)) and R.randint(0,99) < chance: # taper walls
							block = level.blockAt(box.minx+px+x,box.miny+py+y,box.minz+pz+z)
							if block == replaceMaterialID: # Stone (AIR DEBUG)
								setBlock(level,material,box.minx+px+x,box.miny+py+y,box.minz+pz+z)	
		else:
			# Draw a deposit
			size = R.randint(4,8)
			if R.randint(0,100) > 80:
				size = size * 2
			print size
			for y in xrange(-size,size+1):
				for x in xrange(-size,size+1):
					for z in xrange(-size,size+1):
						if y**2+z**2+x**2 <= size**2 and R.randint(0,99) < chance:
							block = level.blockAt(box.minx+px+x,box.miny+py+y,box.minz+pz+z)
							if block == replaceMaterialID: # Stone (AIR DEBUG)
								setBlock(level,material,box.minx+px+x,box.miny+py+y,box.minz+pz+z)	
		
		# Change direction and speed slightly
		theta = theta+R.randint(-3,3)*a*20 # Sharp turns wrap it up tightly
		if R.randint(0,100) > 50:
			phi = phi+R.randint(-1,1)*a*2
		
		if steps % PLUNGE == 0: # MAGICNUMBER: 32 is good. Affects how deep the vertical falls and rises are
			if abs(phi) > 30*a:
				phi = 0+R.randint(-1,1)*a*5 # MAGICNUMBER: How many degrees inclination the tunnel is coming out of a vertical shaft
		velocity = velocity+R.randint(-1,2)
		if velocity < 2:
			velocity = 2
		elif velocity > 3:
			velocity = 3

		#  Move to new segment, repeat
		px = px+velocity*cos(theta)*cos(phi)
		py = py+velocity*sin(phi)
		pz = pz+velocity*sin(theta)*cos(phi)
		
		# Bounds checking
		if px < 16 or pz < 16 or px >= width-16 or pz >= depth-16:
			theta = theta + pi
		if py < 16 or py > height-1:
			phi = phi + pi
		
		steps = steps+1
		if steps > MAXSTEPS:
			keepGoing = False

		# Change the size of the 'brush'
#		if R.randint(0,100) < 20:
#			sizeX = getRandomStep(sizeX,MAXSIZE,R)
#		if R.randint(0,100) < 20:
#			sizeY = getRandomStep(sizeY,2,R)
#		if R.randint(0,100) < 20:
#			sizeZ = getRandomStep(sizeZ,MAXSIZE,R)			
	
	FuncEnd(level,box,options,method)
	
def tunnelTwisting(level,box,options,R,(px,py,pz),MAXLEN,theta,phi,chance):
	method = "Tunnel:Twisting"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start

#	COUNTTUNNELINVOKES = COUNTTUNNELINVOKES + 1
#	if COUNTTUNNELINVOKES > 10: # Too much concurrency will kill you, just as sure as none at all
#		return
	material = (options["Material:"].ID,options["Material:"].blockData)
	PLUNGE = 32
	# replaceMaterialID = options["Block Type to Replace:"].ID
	
	# Direction facing
	if theta == -1:
		theta = a*R.randint(0,MAXANGLES-1)
	if phi == -1:
		phi = a*R.randint(0,MAXANGLES-1)
	velocity = 3
	
	MAXSIZE = 2
	sizeX = 2
	sizeY = 2
	sizeZ = 2

	steps = 0
	MAXSTEPS = MAXLEN
	if MAXLEN < 1:
		MAXSTEPS = R.randint(500,5000) # default for joker input
	keepGoing = True
	while keepGoing == True:
		if R.randint(1,100) < chance:
			# Draw tunnel segment
			for y in xrange(-sizeY,sizeY+1):
				for x in xrange(-sizeX,sizeX+1):
					for z in xrange(-sizeZ,sizeZ+1):
						#if y**2+z**2+x**2 <= 3**2:
						if not ((abs(x) == sizeX and abs(y) == sizeY) or (abs(y) == sizeY and abs(z) == sizeZ) or (abs(x) == sizeX and abs(z) == sizeZ)): # taper walls
							setBlock(level,AIR,box.minx+px+x,box.miny+py+y,box.minz+pz+z)	
		else:
			# Draw a cavern
			size = R.randint(4,8)
			if R.randint(0,100) > 80:
				size = size * 2
			print size
			for y in xrange(-size,size+1):
				for x in xrange(-size,size+1):
					for z in xrange(-size,size+1):
						if y**2+z**2+x**2 <= size**2:
							setBlock(level,AIR,box.minx+px+x,box.miny+py+y,box.minz+pz+z)
			
			# Carve a large cavern
			if R.randint(0,100) > 90: # Increase this likelihood for a big cave
				for k in xrange(0,6):
					tunnelTwisting(level,box,options,R,(px,py,pz),size*2,R.randint(0,MAXANGLES)*a,0,chance)
					
			# Carve a pit bisecting the floor
			if R.randint(0,100) > 10: # Increase this likelihood for more pits
				theta2 = R.random()*pi*2
				for k in xrange(-size,size):
					cosTheta2 = k*cos(theta2)
					sinTheta2 = k*sin(theta2)
					bottom = -R.randint(32,64)
					for y in xrange(bottom,R.randint(0,8*sizeY)):
						for x in xrange(-sizeX,sizeX+1):
							for z in xrange(-sizeZ,sizeZ+1):
								if not ((abs(x) == sizeX and abs(y) == sizeY) or (abs(y) == sizeY and abs(z) == sizeZ) or (abs(x) == sizeX and abs(z) == sizeZ)): # taper walls
									setBlock(level,AIR,box.minx+px+x+cosTheta2,box.miny+py+y,box.minz+pz+z+sinTheta2)	
									if y == bottom and getBlock(level,box.minx+px+x+cosTheta2,box.miny+py+y-1,box.minz+pz+z+sinTheta2) != AIR:
										setBlock(level,(11,0),box.minx+px+x+cosTheta2,box.miny+py+y,box.minz+pz+z+sinTheta2) # LAVA
			# Stalagmites, Tights, Columns, Shawls
			elif R.randint(0,99) > 30:
				theCave = BoundingBox((box.minx+px-size,box.miny+py-size,box.minz+pz-size),(2*size+1,2*size+1,2*size+1))
				columns(level,theCave,options,R,material,R.randint(2,size))
			
		# Change direction and speed slightly
		theta = theta+R.randint(-3,3)*a*20 # Sharp turns wrap it up tightly
		if R.randint(0,100) > 50:
			phi = phi+R.randint(-1,1)*a*2
#		else:
#			phi = phi-R.randint(-1,2)*a*2

		
		if steps % PLUNGE == 0: # MAGICNUMBER: 32 is good. Affects how deep the vertical falls and rises are
			if abs(phi) > 30*a:
				phi = 0+R.randint(-1,1)*a*5 # MAGICNUMBER: How many degrees inclination the tunnel is coming out of a vertical shaft
		velocity = velocity+R.randint(-1,2)
		if velocity < 2:
			velocity = 2
		elif velocity > 3:
			velocity = 3

		# Occasionally branch away
		if steps % 100 == 0: # branch!
			if R.randint(0,100) == 1:
				tunnelTwisting(level,box,options,R,(px,py,pz),int(MAXLEN/2),-1,-1,chance)
		
		#  Move to new segment, repeat
		px = px+velocity*cos(theta)*cos(phi)
		py = py+velocity*sin(phi)
		pz = pz+velocity*sin(theta)*cos(phi)
		
		# Bounds checking
		if px < 16 or pz < 16 or px >= width-16 or pz >= depth-16:
			theta = theta + pi
		if py < 16 or py > height-1:
			phi = phi + pi
		
		steps = steps+1
		if steps > MAXSTEPS:
			keepGoing = False

		# Change the size of the 'brush'
#		if R.randint(0,100) < 20:
#			sizeX = getRandomStep(sizeX,MAXSIZE,R)
#		if R.randint(0,100) < 20:
#			sizeY = getRandomStep(sizeY,2,R)
#		if R.randint(0,100) < 20:
#			sizeZ = getRandomStep(sizeZ,MAXSIZE,R)			
#	COUNTTUNNELINVOKES = COUNTTUNNELINVOKES - 1
	FuncEnd(level,box,options,method)

def tunnelBrownian(level,box,options,R,(px,py,pz),length):
	method = "Tunnel:Brownian"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start
	MAXSIZE = R.randint(2,3)
	
	if length == 0:
		length = R.randint(100,10*R.randint(10,3*(width+depth+height)/3))

	for i in xrange(0,length):
		if R.randint(0,100) < 10:
			(px,py,pz) = (getRandomStep(px,width-1,R),getRandomStep(py,height-1,R),getRandomStep(pz,depth-1,R)) # walk randomly. Or stay still. But... probably walk.
		else:
			(px,py,pz) = (getRandomStep(px,width-1,R),py,getRandomStep(pz,depth-1,R)) # walk randomly. Or stay still. But... probably walk.
		for x in xrange(-1,2):
			for y in xrange(-1,2):
				for z in xrange(-1,2):
					#if not ((abs(x) == sizeX and abs(y) == sizeY) or (abs(y) == sizeY and abs(z) == sizeZ) or (abs(x) == sizeX and abs(z) == sizeZ)): # taper walls
					setBlock(level,AIR,box.minx+px+x,box.miny+py+y,box.minz+pz+z)
		if R.randint(1,1000) == 1:
			tunnelBrownian(level,box,options,R,(px,py,pz),int(length/2))
				
	FuncEnd(level,box,options,method)
	
	return (px,py,pz)
	
def tunnelRandom1(level,box,options,R,(px,py,pz)):
	method = "Tunnel:Random1"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start
	MAXSIZE = R.randint(2,3)
	
	MAXDIR = 7
	
	length = R.randint(100,10*R.randint(10,3*(width+depth+height)/3))
	spread = R.randint(10,40)
	sizeX = R.randint(1,MAXSIZE)
	sizeY = R.randint(1,MAXSIZE)
	sizeZ = R.randint(1,MAXSIZE)
	dir = R.randint(0,MAXDIR+2)
	print length,spread

	for i in xrange(0,length):
		# blit the tunnel shape at the current position, then move off randomly. Repeat
		for x in xrange(-sizeX,sizeX+1):
			for y in xrange(-sizeY,sizeY+1):
				for z in xrange(-sizeZ,sizeZ+1):
					if not ((abs(x) == sizeX and abs(y) == sizeY) or (abs(y) == sizeY and abs(z) == sizeZ) or (abs(x) == sizeX and abs(z) == sizeZ)): # taper walls
						setBlock(level,AIR,box.minx+px+x,box.miny+py+y,box.minz+pz+z)
		if R.randint(0,100) < 20:
			sizeX = getRandomStep(sizeX,MAXSIZE,R)
		if R.randint(0,100) < 20:
			sizeY = getRandomStep(sizeY,2,R)
		if R.randint(0,100) < 20:
			sizeZ = getRandomStep(sizeZ,MAXSIZE,R)
			
		if i%(spread) == 0: # New direction
			if dir <= MAXDIR:
				if R.randint(0,100) > 10:
					dir = dir + R.randint(-1,1)
					if dir < 0:
						dir = MAXDIR
					elif dir > MAXDIR:
						dir = 0
				else:
					dir = R.randint(MAXDIR+1,MAXDIR+2)
			else:
				if R.randint(0,100) > 50:
					dir = R.randint(0,MAXDIR+2)
		
			if dir == 0:
				px = px+1
			elif dir == 1:
				px = px+1
				pz = pz+1
			elif dir == 2:
				pz = pz+1
			elif dir == 3:
				px = px-1
				pz = pz+1
			elif dir == 4:
				px = px-1
			elif dir == 5:
				px = px-1
				pz = pz-1
			elif dir == 6:
				pz = pz-1
			elif dir == 7:
				px = px+1
				pz = pz-1
			elif dir == 8:
				py = py+1
			elif dir == 9:
				py = py-1

				
				
#			if R.randint(0,100) < 2:
#				(px,py,pz) = (getRandomStep(px,width-1,R),getRandomStep(py,height-1,R),getRandomStep(pz,depth-1,R)) # walk randomly. Or stay still. But... probably walk.
#			else:
#				(px,py,pz) = (getRandomStep(px,width-1,R),py,getRandomStep(pz,depth-1,R)) # walk randomly. Or stay still. But... probably walk.

				
	FuncEnd(level,box,options,method)
	
	return (px,py,pz)
	
def getRandomStep(num,MAX,R):
	num = num+R.randint(-1,1)
	if num < 1:
		num = 1
	elif num > MAX:
		num = MAX
	return num

def chaikinSmoothAlgorithm(P): # http://www.idav.ucdavis.edu/education/CAGDNotes/Chaikins-Algorithm/Chaikins-Algorithm.html
	F1 = 0.25
	F2 = 0.75
	Q = []
	(x0,y0,z0) = (-1,-1,-1)
	count = 0
	for (x1,y1,z1) in P:
		if count > 0: # We have a previous point
			(dx,dy,dz) = (x1-x0,y1-y0,z1-z0)
#			Q.append( (x0*F2+x1*F1,0,z0*F2+z1*F1) )
#			Q.append( (x0*F1+x1*F2,0,z0*F1+z1*F2) )

			Q.append( (x0*F2+x1*F1,y0*F2+y1*F1,z0*F2+z1*F1) )
			Q.append( (x0*F1+x1*F2,y0*F1+y1*F2,z0*F1+z1*F2) )

#			Q.append( (x0+dx*F1+*F2,y0*F1+y1*F2,z0*F1+z1*F2) )
#			Q.append( (x0*F2+x1*F1,y0*F2+y1*F1,z0*F2+z1*F1) )
		else:
			count = count+1
		(x0,y0,z0) = (x1,y1,z1)

	return Q
	
def perform(level,box, options):
	''' Feedback to abrightmoore@yahoo.com.au '''
	# Local variables
	method = "Perform"
	(method, (width, height, depth), (centreWidth, centreHeight, centreDepth)) = FuncStart(level,box,options,method) # Log start

	levelCopy = level
	boxCopy = box

	if options["Cache:"] == True:
		levelCopy = level.extractSchematic(box) # Working set
		boxCopy = BoundingBox((0,0,0),(width,height,depth))
	
	underworld(levelCopy,boxCopy,options)

	if options["Cache:"] == True:
		level.copyBlocksFrom(levelCopy, boxCopy, (box.minx, box.miny, box.minz ))
	
	level.markDirtyBox(box)
	FuncEnd(level,box,options,method) # Log end	
	
####################################### LIBS
	
def FuncStart(level, box, options, method):
	# abrightmoore -> shim to prepare a function.
	print '%s: Started at %s' % (method, time.ctime())
	(width, height, depth) = (box.maxx - box.minx, box.maxy - box.miny, box.maxz - box.minz)
	centreWidth = math.floor(width / 2)
	centreHeight = math.floor(height / 2)
	centreDepth = math.floor(depth / 2)	
	# other initialisation methods go here
	return (method, (width, height, depth), (centreWidth, centreHeight, centreDepth))

def FuncEnd(level, box, options, method):
	print '%s: Ended at %s' % (method, time.ctime())
	
def getBoxSize(box):
	return (box.maxx - box.minx, box.maxy - box.miny, box.maxz - box.minz)

def getBlock(level,x,y,z):
	return (level.blockAt(int(x),int(y),int(z)), level.blockDataAt(int(x),int(y),int(z)))

def setBlock(level, (block, data), x, y, z):
	level.setBlockAt(int(x), int(y), int(z), block)
	level.setBlockDataAt(int(x), int(y), int(z), data)

def drawTriangleEdge(level, box, options, (p1x, p1y, p1z), (p2x, p2y, p2z), (p3x, p3y, p3z), materialEdge):
	drawLine(level, materialEdge, (p1x, p1y, p1z), (p2x, p2y, p2z) )
	drawLine(level, materialEdge, (p1x, p1y, p1z), (p3x, p3y, p3z) )
	drawLine(level, materialEdge, (p2x, p2y, p2z), (p3x, p3y, p3z) )
	
# Ye Olde GFX Libraries
def drawLine(scratchpad, (blockID, blockData), (x,y,z), (x1,y1,z1) ):
	drawLineConstrained(scratchpad, (blockID, blockData), (x,y,z), (x1,y1,z1), 0 )

def drawLine1(scratchpad, (blockID, blockData), (x,y,z), (x1,y1,z1) ):
	for px, py, pz in bresenham.bresenham((x,y,z),(x1,y1,z1)):
		setBlock(scratchpad,(blockID, blockData),px,py,pz)
	setBlock(scratchpad,(blockID, blockData),x1,y1,z1)
	
def drawLineConstrained(scratchpad, (blockID, blockData), (x,y,z), (x1,y1,z1), maxLength ):
	dx = x1 - x
	dy = y1 - y
	dz = z1 - z

	distHoriz = dx*dx + dz*dz
	distance = sqrt(dy*dy + distHoriz)

	if distance < maxLength or maxLength < 1:
		phi = atan2(dy, sqrt(distHoriz))
		theta = atan2(dz, dx)

		iter = 0
		while iter <= distance:
			scratchpad.setBlockAt((int)(x+iter*cos(theta)*cos(phi)), (int)(y+iter*sin(phi)), (int)(z+iter*sin(theta)*cos(phi)), blockID)
			scratchpad.setBlockDataAt((int)(x+iter*cos(theta)*cos(phi)), (int)(y+iter*sin(phi)), (int)(z+iter*sin(theta)*cos(phi)), blockData)
			iter = iter+0.5 # slightly oversample because I lack faith.

def drawLineConstrainedRandom(scratchpad, (blockID, blockData), (x,y,z), (x1,y1,z1), frequency ):
	dx = x1 - x
	dy = y1 - y
	dz = z1 - z

	distHoriz = dx*dx + dz*dz
	distance = sqrt(dy*dy + distHoriz)


	phi = atan2(dy, sqrt(distHoriz))
	theta = atan2(dz, dx)

	iter = 0
	while iter <= distance:
		if randint(0,99) < frequency:
			scratchpad.setBlockAt((int)(x+iter*cos(theta)*cos(phi)), (int)(y+iter*sin(phi)), (int)(z+iter*sin(theta)*cos(phi)), blockID)
			scratchpad.setBlockDataAt((int)(x+iter*cos(theta)*cos(phi)), (int)(y+iter*sin(phi)), (int)(z+iter*sin(theta)*cos(phi)), blockData)
		iter = iter+0.5 # slightly oversample because I lack faith.

def drawTriangle(level, (p1x, p1y, p1z), (p2x, p2y, p2z), (p3x, p3y, p3z), materialEdge, materialFill):
	if materialFill != (0,0):
		# for each step along the 'base' draw a line from the apex
		dx = p3x - p2x
		dy = p3y - p2y
		dz = p3z - p2z

		distHoriz = dx*dx + dz*dz
		distance = sqrt(dy*dy + distHoriz)
		
		phi = atan2(dy, sqrt(distHoriz))
		theta = atan2(dz, dx)

		iter = 0
		while iter <= distance:
			(px, py, pz) = ((int)(p2x+iter*cos(theta)*cos(phi)), (int)(p2y+iter*sin(phi)), (int)(p2z+iter*sin(theta)*cos(phi)))
			
			iter = iter+0.5 # slightly oversample because I lack faith.
			drawLine(level, materialFill, (px, py, pz), (p1x, p1y, p1z) )
	
	
	drawLine(level, materialEdge, (p1x, p1y, p1z), (p2x, p2y, p2z) )
	drawLine(level, materialEdge, (p1x, p1y, p1z), (p3x, p3y, p3z) )
	drawLine(level, materialEdge, (p2x, p2y, p2z), (p3x, p3y, p3z) )
