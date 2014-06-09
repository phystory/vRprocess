import os
import numpy as np
import struct
import binascii
import array
from matplotlib  import pyplot as plt
from matplotlib import animation as animation

"""
Python code to read tso file based on matlab code
function [zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx] = read_ts_file( ts_filename )
%--------------------------------------------------------------------------
%[zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx] = read_ts_file( ts_filename ) 
% Reads XNet binary output file.
% Inputs>  ts_filename: name of binary file 
% Outputs< zz: proton number for each isotope in the network.
%          aa: total nucleon number for each isotope in the network
%          xmf: time series of mass fractions 
%          time: discrete temporal evolution
%          temperature:  time series of temperature
%          density: time series of density
%          timestep:  time series of discrete timesteps 
%          edot: time series of energy generation 
%          flux_end, the starting and ending points for each reaction flux.
%          flux: time seris of the reaction fluxes in the network.
%--------------------------------------------------------------------------
"""

length = os.path.getsize("tso1")
desc1 = []
desc2 = []
desc3 = []
data_desc = []
abund_file = []
abund_desc = []
thermo_file = []
thermo_desc = []
zz = []
aa = []
xmf = []
flx = [[3]]
flx_end = []
time = []
temperature = []
density = []
timestep = []
edot = []
with open("tso1",'rb') as file_id:
	record_length1=struct.unpack('i',file_id.read(4))[0]
	for ith in range (0,80):
		desc1.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		desc2.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		desc3.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		data_desc.append(struct.unpack('c',file_id.read(1))[0])
	record_length2=struct.unpack('i',file_id.read(4))[0]
	record_length1=struct.unpack('i',file_id.read(4))[0]
	kstmx    =struct.unpack('i',file_id.read(4))[0]
	kitmx    =struct.unpack('i',file_id.read(4))[0]
	iweak    =struct.unpack('i',file_id.read(4))[0]
	iscrn    =struct.unpack('i',file_id.read(4))[0]
	iconvc   =struct.unpack('i',file_id.read(4))[0]
	changemx =struct.unpack('@d',file_id.read(8))[0]
	tolm     =struct.unpack('@d',file_id.read(8))[0]
	tolc     =struct.unpack('@d',file_id.read(8))[0]
	yacc     =struct.unpack('@d',file_id.read(8))[0]
	ymin     =struct.unpack('@d',file_id.read(8))[0]
	tdel_mm  =struct.unpack('@d',file_id.read(8))[0]
	record_length2=struct.unpack('i',file_id.read(4))[0]

	# Read Abundance Info
	record_length1=struct.unpack('i',file_id.read(4))[0]
	for ith in range (0,80):
		abund_file.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		abund_desc.append(struct.unpack('c',file_id.read(1))[0])
	record_length2=struct.unpack('i',file_id.read(4))

	# Read Thermodynamic Info
	record_length1=struct.unpack('i',file_id.read(4))[0]
	for ith in range (0,80):
		thermo_file.append(struct.unpack('c',file_id.read(1))[0])
	for ith in range (0,80):
		thermo_desc.append(struct.unpack('c',file_id.read(1))[0])
	record_length2=struct.unpack('i',file_id.read(4))[0]

	# Read Nuclear Info
	record_length1=struct.unpack('i',file_id.read(4))[0]
	ny            =struct.unpack('i',file_id.read(4))[0]
	print "ny=",ny
	for ith in range (0,ny):
		zz.append(struct.unpack('@d',file_id.read(8))[0])
	for ith in range (0,ny):
		aa.append(struct.unpack('@d',file_id.read(8))[0])
	record_length2=struct.unpack('i',file_id.read(4))[0]
	#print zz,aa

	# Read Flux Info
	record_length1=struct.unpack('i',file_id.read(4))[0]
	nflx          =struct.unpack('i',file_id.read(4))[0]
	if nflx > 0:
		for ith in range(0,nflx):
			flx_end1=struct.unpack('i',file_id.read(4))[0]
			flx_end2=struct.unpack('i',file_id.read(4))[0]
			#print "flx_end1=",flx_end1,"flx_end2=",flx_end2
			#flx_end[ith].append(struct.unpack('i',file_id.read(4))[0])
			#flx_end[ith].append(struct.unpack('i',file_id.read(4))[0])
	record_length2=struct.unpack('i',file_id.read(4))[0]

	# Read data from each timestep
	for k in range(1, kstmx):
		record_length1 =struct.unpack('i',file_id.read(4))[0]
		# If end of file, exit 
		if not record_length1:
			break
		# Otherwise read data
		kstep          =struct.unpack('i',file_id.read(4))[0]
		#time.append(struct.unpack('d',file_id.read(8))[0])
		#temperature.append(struct.unpack('d',file_id.read(8))[0])
		#density.append(struct.unpack('d',file_id.read(8))[0])
		#timestep.append(struct.unpack('d',file_id.read(8))[0])
		#edot.append(struct.unpack('d',file_id.read(8))[0])
		timel = struct.unpack('@d',file_id.read(8))[0]
		temperaturel = struct.unpack('d',file_id.read(8))[0]
		densityl = struct.unpack('@d',file_id.read(8))[0]
		timestepl = struct.unpack('@d',file_id.read(8))[0]
		edotl = struct.unpack('@d',file_id.read(8))[0]
		for ith in range(0,ny):
			xmfl = struct.unpack('@d',file_id.read(8))[0]
			print k,timel,aa[ith],zz[ith],xmfl
		if nflx > 0: 
			for ith in range(0,nflx):
				flxl = struct.unpack('@d',file_id.read(8))[0]
				#print "flxl=",flxl
		record_length2 =struct.unpack('i',file_id.read(4))[0]
"""
		for ith in range(0,ny):
			xmf[ith].append([])
			xmf[ith].append(struct.unpack('d',file_id.read(8))[0])
		if nflx > 0: 
			for ith in range(0,ny):
				flx[ith].append([])
				flx[ith].append(struct.unpack('i',file_id.read(8))[0])
"""
	#print kstep,time,temperature,density
	# Convert abundance to mass fraction
	#for n in range(1, ny):
		#xmf[n] = aa[n]*xmf[n]

	# Calculate # of timesteps
	#nstep = len(time,2)

"""
% Read Run Descriptions
  record_length1=fread(file_id,1,'int32');
  desc1    =setstr(fread(file_id,80,'uchar'));
  desc2    =setstr(fread(file_id,80,'uchar'));
  desc3    =setstr(fread(file_id,80,'uchar'));
  data_desc=setstr(fread(file_id,80,'uchar'));
  record_length2=fread(file_id,1,'int32');

% Read Run Settings
  record_length1=fread(file_id,1,'int32');
  kstmx    =fread(file_id,1,'int32');
  kitmx    =fread(file_id,1,'int32');
  iweak    =fread(file_id,1,'int32');
  iscrn    =fread(file_id,1,'int32');
  iconvc   =fread(file_id,1,'int32');
  changemx =fread(file_id,1,'float64');
  tolm     =fread(file_id,1,'float64');
  tolc     =fread(file_id,1,'float64');
  yacc     =fread(file_id,1,'float64');
  ymin     =fread(file_id,1,'float64');
  tdel_mm  =fread(file_id,1,'float64');
  record_length2=fread(file_id,1,'int32');

% Read Abundance Info
  record_length1=fread(file_id,1,'int32');
  abund_file =setstr(fread(file_id,80,'uchar'));
  abund_desc =setstr(fread(file_id,80,'uchar'));
  record_length2=fread(file_id,1,'int32');

% Read Thermodynamic Info
  record_length1=fread(file_id,1,'int32');
  thermo_file =setstr(fread(file_id,80,'uchar'));
  thermo_desc =setstr(fread(file_id,80,'uchar'));
  record_length2=fread(file_id,1,'int32');

% Read Nuclear Info
  record_length1=fread(file_id,1,'int32');
  ny            =fread(file_id,1,'int32')
  zz            =fread(file_id,ny,'float64');
  aa            =fread(file_id,ny,'float64');
  record_length2=fread(file_id,1,'int32');

% Read Flux Info
  record_length1=fread(file_id,1,'int32');
  nflx            =fread(file_id,1,'int32')
  if nflx>0
    flx_end(:,1)    =fread(file_id,nflx,'int32');
    flx_end(:,2)    =fread(file_id,nflx,'int32');
  end
  record_length2=fread(file_id,1,'int32');

% Read data from each timestep
  for k= 1:kstmx
    record_length1 =fread(file_id,1,'int32');
  % If end of file, exit 
    if isempty(record_length1)
      break
    end
  % Otherwise read data
    kstep          =fread(file_id,1,'int32');
    time(k)        =fread(file_id,1,'float64');
    temperature(k) =fread(file_id,1,'float64');
    density(k)     =fread(file_id,1,'float64');
    timestep(k)    =fread(file_id,1,'float64');
    edot(k)        =fread(file_id,1,'float64');
    xmf(:,k)       =fread(file_id,ny,'float64');
    if nflx>0 
      flx(:,k)       =fread(file_id,nflx,'float64');
    end
    record_length2 =fread(file_id,1,'int32');
  end

% Convert abundance to mass fraction
  for n=1:ny
    xmf(n,:) = aa(n).*xmf(n,:);
  end

% Calculate # of timesteps
  nstep = size(time,2)
  
end
"""

