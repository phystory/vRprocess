#!/usr/bin/python

#Auslesen von Daten-Files
def read_ascii(file,columns,lines=None,dtype=float):	     	 
    """[X,Y,Z] = READ('FILE',[1,4,13],lines=[10,1000])	     	 
    Read columns 1,4 and 13 from 'FILE'  from line 10 to 1000	 
    into array X,Y and Z				     	 
    """ 						     	 
    import sys,os,string,re
    import numpy as N
    try: import mmap					     	 
    except: pass					     	 
    f = open(file,'r+') 				     	 
    f.seek(0,2) 					     	 
    fsize = f.tell()					     	 
    f.seek(0,0) 					     	 
    start = 1						     	 
    end = 0						     	 
    numcols=len(columns)				     	    							     	 
    if lines is None:					     	 
	lines = [start,end]				     	 
    elif len(lines) != 2:				     	 
	raise "lines!=[start,end]"			     	 
    [start,end] = map(int,lines)			     	 
    start = start - 1					     	 							     	 
    wc = string.split(os.popen("wc "+file).read())	     	 
    if (end == 0) or (end > int(wc[0])):		     	 
	end = int(wc[0])				     	 
    numlines = end - start				     	     							     	 
    block = 65536  # chunk to read, in bytes		     	 
    data = mmap.mmap(f.fileno(), fsize) 		     	 
    for i in range(start):				     	 
	junk = data.readline()				     	 							     	 
    numcols =  len(columns)
#    print numcols,numlines				     	 
    if dtype == int:					     	 
       a = N.zeros((numlines,numcols), 'l')	     	 
    else:						     	 
       a = N.zeros((numlines,numcols), 'd')	     	 
    i = 0						     	 
    line = data.readline()				     	 
    while line and ( i < numlines):			     	 
	line=re.sub('[+,-][0-9][0-9][0-9] ','e-99 ',line)
	line=re.sub('\*','0',line)
	d = N.array(map(dtype,string.split(line)))    
	a[i,:] = N.take(d,columns)		     	 
	line = data.readline()				     	 
	i = i + 1					     	 
	# writing reading status 
        percent = int(i*100/numlines)
        sys.stdout.write("\rreading file " + "...%d%%" % percent)
        sys.stdout.flush()
    data.close()					     	 
    f.close()						     	 
    b = []						     	 
    print ''
    for i in range(numcols):				     	 
	b.append(a[:,i])				     	 
    del a						     	 
    return  b	 

