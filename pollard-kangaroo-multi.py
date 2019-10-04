#!/usr/bin/python

# based on code by 57fe, 2019
# fe57.org/forum/thread.php?board=4&thema=1#1

#######################
# print() compatibility python 2/3
from __future__ import print_function
#######################
# users settings

pow2bits	= 42	# bits (suborder) range search of keyspace (expected location of privkey)


Ntimeit		= 10		# times for avg runtime
timeit_eachnewprvkey = True	# gen new privkey each loop?


flag_profile	= "byPollard"		# best, expected 2w^(1/2)/cores jumps
#flag_profile	= "NaiveSplitRange"	# norm, expected 2(w/cores)^(1/2) jumps
#flag_profile	= "byOorschot&Wiener"	# possibly bad implementation, low efficiency
#flag_profile	= "stupid_random"	# test for compare, not recommended

#######################
# service settings

max_cpu_cores	= 128
min_cpu_cores	= 1

flag_verbose	= 0	# 0, 1, 2

prngseed	= 0	# 0 for random, or any for replay results

version = '1.04'

# low order pubkeys
# default_table (demo/debug)
pubkeys = {
	  16: ('029d8c5d35231d75eb87fd2c5f05f65281ed9573dc41853288c62ee94eb2590b7a', 0xc936)
	, 24: ('036ea839d22847ee1dce3bfc5b11f6cf785b0682db58c35b63d1342eb221c3490c', 0xdc2a04)
	, 32: ('0209c58240e50e3ba3f833c82655e8725c037a2294e14cf5d73a5df8d56159de69', 0xb862a62e)
	, 33: ('02ed949eaca31df5e8be9bf46adc1dfae1734b8900dcc303606831372955c728da', False) #0x01abcd1234
	, 40: ('03a2efa402fd5268400c77c20e574ba86409ededee7c4020e4b9f0edbee53de0d4', 0xe9ae4933d6)
	, 50: ('03f46f41027bbf44fafd6b059091b900dad41e6845b2241dc3254c7cdd3c5a16c6', 0x022bd43c2e9354)
	, 55: ('0385a30d8413af4f8f9e6312400f2d194fe14f02e719b24c3f83bf1fd233a8f963', 0x6abe1f9b67e114)
	, 60: ('0348e843dc5b1bd246e6309b4924b81543d02b16c8083df973a89ce2c7eb89a10d', 0x0FC07A1825367BBE)
	, 70: ('0290e6900a58d33393bc1097b5aed31f2e4e7cbd3e5466af958665bc0121248483', 0x349B84B6431A6C4EF1)
	, 80: ('037e1238f7b1ce757df94faa9a2eb261bf0aeb9f84dbf81212104e78931c2a19dc', 0xEA1A5C66DCC11B5AD180)
	, 90: ('035c38bd9ae4b10e8a250857006f3cfd98ab15a6196d9f4dfd25bc7ecc77d788d5', 0x02CE00BB2136A445C71E85BF)
	,100: ('03d2063d40402f030d4cc71331468827aa41a8a09bd6fd801ba77fb64f8e67e617', 0x0af55fc59c335c8ec67ed24826)
	,105: ('03bcf7ce887ffca5e62c9cabbdb7ffa71dc183c52c04ff4ee5ee82e0c55c39d77b', False)
}

#######################
# import

import os
import sys
import time
import math
import random

# gmpy2 is the fastest!
# download file .whl from https://www.lfd.uci.edu/~gohlke/pythonlibs/
# [windows>python.exe - m] pip install gmpy2-2.0.8-cp37-cp37m-win_amd64.whl
try:
	# https://www.lfd.uci.edu/~gohlke/pythonlibs/
	import gmpy2
except:
	flag_gmpy2 = False
	print("[warn] lib gmpy2 not found. full speed is not achievable!")
else:
	flag_gmpy2 = True


# debug lib
if 0:
	flag_gmpy2 = 0


import multiprocessing as mp


#######################
# python2/3 compatibility

#import sys
#import time
if sys.version_info[0] == 2:
	from time import clock
else:
	from time import perf_counter
	from time import process_time
	clock = time.perf_counter
	xrange=range
	raw_input=input


#######################
# ec secp256k1

A_short	= 0
B_short	= 7
#modulo	= 2**256-2**32-2**9-2**8-2**7-2**6-2**4-1
modulo	= 115792089237316195423570985008687907853269984665640564039457584007908834671663
order	= 115792089237316195423570985008687907852837564279074904382605163141518161494337
#modulo	= 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
#order	= 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
Gx	= 55066263022277343669578718895168534326250603453777594175500187360389116729240
Gy	= 32670510020758816978083085130507043184471273380659243275938904335757337482424
#Gx	= 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
#Gy	= 0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8


# python2+gmpy2 speed-up +8%
if flag_gmpy2:
	A_short	= gmpy2.mpz(A_short)
	B_short	= gmpy2.mpz(B_short)
	modulo	= gmpy2.mpz(modulo)
	order	= gmpy2.mpz(order)
	Gx	= gmpy2.mpz(Gx)
	Gy	= gmpy2.mpz(Gy)


class Point:
	def __init__(self, x=0, y=0):		# Affine
	#def __init__(self, x=0, y=0, z=1):	# Jacobian
		self.x = x
		self.y = y
		#self.z = 1			# Jacobian

Gp = Point(Gx,Gy) 
Zp = Point(0,0)	# zero-point, infinite in real x,y - plane


#######################
# math, raw python


# from arulberoEC library
# more fastest
def invert(b, p=modulo):	
	u, v = b%p, p
	x1, x2 = 1, 0
	while u != 1:
		#q = v//u
		#r = v-q*u
		q, r = divmod(v,u)		
		x = x2-q*x1
		v = u
		u = r
		x2 = x1
		x1 = x
	return x1%p


#######################
# Affine coordinates (X,Y,Z=1)

# specific of python: x*x... more faster than x**2 !!!
# ..so option "0S" (without squaring) more preferred


# A + A -> A (1I, 2M, 1S)
# A + A -> A (1I, 3M, 0S)
def add_a(A, B, p=modulo):
	R = Point()
	dx = B.x - A.x
	dy = B.y - A.y	
	if flag_gmpy2:
		c = dy * gmpy2.invert(dx, p) % p	# 1I,1M
	else:
		c = dy * invert(dx, p) % p		# 1I,1M

	#R.x = (c**2 - A.x - B.x) % p	# 0M,1S
	R.x = (c*c - A.x - B.x) % p	# 1M,0S
	#R.x = (int(math.pow(c,2)) - A.x - B.x) % p	# slow

	R.y = (c*(A.x - R.x) - A.y) % p	# 1M
	return R


# 2 * A -> A (1I, 5M, 2S)
# 2 * A -> A (1I, 7M, 0S)
# 2 * A -> A (1I, 4M, 0S)
def mul_2a(A, p=modulo):
	R = Point()
	if flag_gmpy2:
		#c = 3 * A.x**2 * gmpy2.invert(2*A.y, p) % p	# 1I,3M,1S
		#c = 3 * A.x * A.x * gmpy2.invert(2*A.y, p) % p	# 1I,4M,0S

		c = A.x * A.x * gmpy2.invert(A.y+A.y, p) 	# 1I,2M,0S
		c = (c + c + c) % p;
	else:
		#c = 3 * A.x**2 * invert(2*A.y, p) % p		# 1I,3M,1S
		#c = 3 * A.x * A.x * invert(2*A.y, p) % p	# 1I,4M,0S

		c = A.x * A.x * invert(A.y+A.y, p)		# 1I,2M,0S
		c = (c + c + c) % p;

	#R.x = (c**2 - 2*A.x) % p	# 1M,1S
	#R.x = (c*c - 2*A.x) % p	# 2M,0S
	R.x = (c*c - A.x - A.x) % p	# 1M,0S

	R.y = (c*(A.x - R.x) - A.y) % p	# 1M
	return R


# k * A -> A
def mul_ka(k, A=Gp, p=modulo):
	if k == 0: return Zp
	elif k == 1: return A
	elif (k%2 == 0):
		return mul_ka(k//2, mul_2a(A, p), p)
	else:
		return add_a(A, mul_ka( (k-1)//2, mul_2a(A, p), p), p)


#######################
# support functions

# calculation Y from X if pubkey is compressed
# more fastest
def getX2Y(X, y_parity, p=modulo):
	
	Y = 3
	tmp = 1
	while Y:
		if Y & 1:
			tmp = tmp*X % p
		Y >>= 1
		X = X*X % p

	X = (tmp+7) % p

	Y = (p+1)//4
	tmp = 1
	while Y:
		if Y & 1:
			tmp = tmp*X % p
		Y >>= 1
		X = X*X % p

	Y = tmp

	if Y%2 != y_parity:
		Y = -Y % p

	return Y


def save2file(path, mode, data):
	fp = open(path, mode)
	if type(data) in (list,tuple,dict,set):
		fp.writelines(data)
	else:
	#elif type(data) in (str,int):
		fp.write(data)
	fp.close()


def usage(bits=32):
	print('[usage] %s [bits] [pubkey]'%(sys.argv[0]))
	print('        %s %s'%(sys.argv[0],bits))
	print('        %s %s %s'%(sys.argv[0],bits,pubkeys[bits][0]))
	print('        %s 12ABCDEF:FFFF0000 %s'%(sys.argv[0],pubkeys[bits][0]))
	exit(-1)


def prefSI(num):
	prefSI_index = 0
	# Kilo/Mega/Giga/Tera/Peta/Exa/Zetta/Yotta
	dict_prefSI = {0:'', 1:'K', 2:'M', 3:'G', 4:'T', 5:'P', 6:'E', 7:'Z', 8:'Y'}
	num *= 1.0
	while( int(num/1000) > 0): 
		prefSI_index += 1
		num /= 1000
	if prefSI_index >= len(dict_prefSI):
		return ('infini')
	else:
		return ('%.1f'%num)+dict_prefSI[prefSI_index]
#print('%s' % prefSI(int(sys.argv[1])));exit(1)


def time_format(time, v=(1,1,1,1,1,1,0,0)):
	sec  = int(time)
	msec = int((time%1)*1000)
	usec = int((((time%1)*1000)%1)*1000)
	res  = ''	
	if v[0]: 
		Y_tmp	= (sec//(60*60*24*30))//12
		if Y_tmp>0: 
			res += ' '+'%s'%(Y_tmp if Y_tmp<10**3 else prefSI(Y_tmp))		+'y'	# year
	if v[1]:
		M_tmp	= (sec//(60*60*24*30))%12
		if M_tmp>0 or Y_tmp>0: 
			res += ' '+'%02s'%str(M_tmp)						+'m'	# month
	if v[2]: 
		d_tmp	= (sec//(60*60*24))%30
		if M_tmp>0 or Y_tmp>0 or d_tmp:
			res += ' '+'%02s'%str(d_tmp)						+'d'	# day
	if v[3]:
		h_tmp	= (sec//(60*60))%24
		if 1:
			res += ' '+'%02d'%int(h_tmp)						+''	# hour
	if v[4]: 
		m_tmp	= (sec//(60*1))%60
		if 1:
			res += ':'+'%02d'%int(m_tmp)						+''	# min
	if v[5]: 
		s_tmp	= (sec//(1*1))%60
		if 1:
			res += ':'+'%02d'%int(s_tmp)						+'s'	# sec
	if v[6]: 
		ms_tmp	= msec
		if v[6]==1:
			res += ' '+'%03d'%(ms_tmp)
		elif Y_tmp==M_tmp==d_tmp==h_tmp==m_tmp==s_tmp==0: 
			res += ' '+'%03d'%(ms_tmp)						+'ms'	# msec
	if v[7]:
		us_tmp	= usec
		if v[6]==1:
			res += ' '+'%03d'%(us_tmp)
		elif Y_tmp==M_tmp==d_tmp==h_tmp==m_tmp==s_tmp==ms_tmp==0: 
			res += ' '+'%03d'%(us_tmp)						+'us'	# usec
	return res
#print('[time] %s'%time_format(int(sys.argv[1]), (1,1,1,1,1,1,0,0)));exit(1)


# # 1<<123 === 2**123, its same, byte shift trick, but 1<< is more x10 faster!
def benchmark_pow2(pow2max=9999):
	tmp=0
	t0 = time.time()
	for i in xrange(1,pow2max):
		tmp += 1<<i
	time1 = time.time()-t0
	print('[%s] %ssec' % ('1<<', time1))

	tmp=0
	t0 = time.time()
	for i in xrange(1,pow2max):
		tmp += 2**i
	time2 = time.time()-t0
	print('[%s] %ssec' % ('2**', time2))

	print('[1<<] %.0f faster than [2**]' % (time2/time1) )
#benchmark_pow2();exit(1)


# fast get X coordinate from point
def getXcoord(itpoint):
	if flag_gmpy2:
		Xcoord = itpoint.x
		#Ycoord = itpoint.y
	else:
		Xcoord = itpoint.x
		#Ycoord = itpoint.y
	return Xcoord #,Ycoord


# get hex pubkey from int prvkey
def getPubkey(new_prvkey, flag_compress):
	if flag_gmpy2:
		Ptmp = mul_ka(new_prvkey)
		Xcoord = Ptmp.x
		Ycoord = Ptmp.y
	else:
		Ptmp = mul_ka(new_prvkey)
		Xcoord = Ptmp.x
		Ycoord = Ptmp.y

	if flag_compress:
		if (Ycoord % 2) == 0:
			new_pubkey = '02%064x' % int(hex(Xcoord)[2:66],16)
		else:
			new_pubkey = '03%064x' % int(hex(Xcoord)[2:66],16)
	else:
		new_pubkey = '04%064x%064x' % (int(hex(Xcoord)[2:66],16), int(hex(Ycoord)[2:66],16))

	return new_pubkey


# get JmaxofSp
def getJmaxofSp(optimalmeanjumpsize, dS):
	if flag_verbose > 0: 
		print('[optimal_mean_jumpsize] %s' % optimalmeanjumpsize)

	sumjumpsize = 0

	for i in range(1,len(dS)):

		#sumjumpsize = (2**i)-1
		#sumjumpsize += 2**(i-1)
		sumjumpsize += dS[i-1]

		now_meanjumpsize	= int(round(1.0*(sumjumpsize)/(i)))

		#next_meanjumpsize	= int(round(1.0*(sumjumpsize+2**i)/(i+1)))
		next_meanjumpsize	= int(round(1.0*(sumjumpsize+dS[i])/(i+1)))

		if flag_verbose > 1: 
			print('[meanjumpsize#Sp[%d]] %s(now) <= %s(optimal) <= %s(next)' % (i, now_meanjumpsize, optimalmeanjumpsize, next_meanjumpsize ))


		if  optimalmeanjumpsize - now_meanjumpsize <= next_meanjumpsize - optimalmeanjumpsize : 
			if flag_verbose > 0: 
				print('[meanjumpsize#Sp[%d]] %s(now) <= %s(optimal) <= %s(next)' % (i, now_meanjumpsize, optimalmeanjumpsize, next_meanjumpsize ))

			# location in keyspace on the strip
			if flag_verbose > 0:
				if (optimalmeanjumpsize - now_meanjumpsize) >= 0:
					len100perc = 60
					size1perc = (next_meanjumpsize-now_meanjumpsize)//len100perc
					print("[i] Sp[%s]|%s%s%s|Sp[%s]" % ( i
						, '-'*(abs(optimalmeanjumpsize - now_meanjumpsize)//size1perc)
						, 'J'
						, '-'*(abs(next_meanjumpsize - optimalmeanjumpsize)//size1perc)
						, i+1 )
					)
					if 1.0*abs(optimalmeanjumpsize-now_meanjumpsize)/abs(next_meanjumpsize-optimalmeanjumpsize) >= 0.25 :
						print("[i] this Sp set has low efficiency (over -25%) for this mean jumpsize")
				else:
					# recovery last step
					now_meanjumpsize	= int(round(1.0*(sumjumpsize-dS[i-1])/(i-1)))
					next_meanjumpsize	= int(round(1.0*(sumjumpsize)/(i)))

					len100perc = 60
					size1perc = (next_meanjumpsize-now_meanjumpsize)//len100perc
					print("[i] Sp[%s]|%s%s%s|Sp[%s]" % ( i-1
						, '-'*(abs(optimalmeanjumpsize - now_meanjumpsize)//size1perc)
						, 'J'
						, '-'*(abs(next_meanjumpsize - optimalmeanjumpsize)//size1perc)
						, i )
					)
					if 1.0*abs(next_meanjumpsize-optimalmeanjumpsize)/abs(optimalmeanjumpsize-now_meanjumpsize) >= 0.25 :
						print("[i] this Sp set has low efficiency (over -25%) for this mean jumpsize")
				#exit(1)

			if flag_verbose > 0: 
				print('[JmaxofSp] Sp[%s]=%s nearer to optimal mean jumpsize of Sp set' % (i, now_meanjumpsize))

			return i

	print("\n[FATAL_ERROR] JmaxofSp not defined!\n"); exit(-1)


# Checks whether the given point lies on the elliptic curve
def is_on_curve(Xcoord,Ycoord, p=modulo):
	# convert short->full:  A_full=0, B_full=A_short, C_full=B_short
	# convert full->short:  A_short=B_full, B_short=C_full (if A_full!=0 - convert impossible!)

	# short form Weierstrass cubic 
	# A_short, B_short
	# y^2 = x^3 + a*x + b over Fp
        return ((Ycoord * Ycoord) - (Xcoord * Xcoord * Xcoord) - (A_short * Xcoord) - B_short) % p == 0

	# full  form Weierstrass cubic 
	# A_full, B_full, C_full
	# y^2 = x^3 + a*x^2 + b*x + c over Fp
        #return ((Ycoord * Ycoord) - (Xcoord * Xcoord * Xcoord) - (A_full * Xcoord * Xcoord) - (B_full * Xcoord) - C_full) % p == 0


#######################
# KANGAROO

def KANGAROO(id_uniq, Sp, dS, Kp, dK, DPmodule, JmaxofSp, jump_step, send2parent, recv_repair):

	#parrent_pid = os.getppid()
	child_pid  = os.getpid()
	child_name = mp.current_process().name
	if flag_verbose > 0: 
		print("[childs][%s#%s] run.." % (child_name, child_pid))


	flag_delay = "by time"
	#flag_delay = "by jumps"


	time_step = 0.2
	t0 = t1 = t2 = time.time()

	prvkey = False;
	countj = last_countj = 0
	repair_Kp, repair_dK = list(), list()

	# main loop
	while (1):
		
		# kangaroo
		if 1:
			countj += 1

			# Xcoord
			Xcoord = getXcoord(Kp)

			pw = Xcoord % JmaxofSp
			pw = int(pw)
			#nowjumpsize = 1<<pw
			nowjumpsize = dS[pw]

			# check, is it distinguished point?
			if Xcoord % DPmodule == 0:

				# send new distinguished point to parent
				send2parent.put_nowait({'id_uniq': id_uniq, 'pid': child_pid, 'name': child_name
						, 'diffjumps': False, 'Xcoord': Xcoord, 'dK': dK})
				repair_Kp.append(Kp)
				repair_dK.append(dK)

			dK += nowjumpsize
			Kp = add_a(Kp, Sp[pw])


		# stat and repair
		if not (countj % jump_step) and countj != 0:

			t2 = time.time()

			# stat
			if flag_delay == "by jumps" or (t2-t1)>=time_step:

				# send diffjumps to parent
				send2parent.put_nowait(
						{'id_uniq': id_uniq, 'pid': child_pid, 'name': child_name
						, 'diffjumps': countj-last_countj, 'Xcoord': False, 'dK': False}
				)

				t1 = t2
				last_countj = countj

				# repair
			if not recv_repair.empty():
				#if 0:
					# recv msg from parent about repair 
					msg = recv_repair.get_nowait()
					dK		= msg['dK']

					# check if already repaired
					try:
						index = repair_dK.index(dK)
					except:
						pass
					else:
						# generate random even offset
						pow2offset = random.randint(0, JmaxofSp)
						if pow2offset == 0:
							dK = repair_dK[index] + 1
							Kp = add_a(repair_Kp[index], Sp[0])
						else:
							dK = repair_dK[index] + (1<<pow2offset) +1
							Kp = add_a(add_a(repair_Kp[index], Sp[pow2offset]), Sp[0])

						if flag_verbose > 0: 
							printstr  = '\n'
							printstr += '[childs][%s#%s] repair#*/*:' % (child_name, child_pid)
							#printstr += ' 0x%064x' % (getXcoord(repair_Kp[index]))
							printstr += ' 0x%x' % (repair_dK[index])
							#printstr += '%65s' % (' ')
							print(printstr)


						#repair_Kp.clear();repair_dK.clear();
						repair_Kp, repair_dK = list(), list()
						


#######################
#main

if __name__ == '__main__':

	#print('[os] %s' % os.name)
	if os.name == 'nt':
		mp.freeze_support()

	try:
		mp.set_start_method('spawn')
	except:
		pass

                                       ##
	print("[################################################]")
	print("[#    Pollard-kangaroo PrivKey Recovery Tool    #]")
	print("[#           ecdsa on curve secp256k1           #]")
	print("[#                   multicore                  #]");
	print("[#                    ver%04s                   #]"%version);
	print("[################################################]")

	if len(sys.argv) > 1 and str(sys.argv[1]) in ('--help','-h','/?') :
		usage()

	print('[date] %s' % (time.ctime()))
	print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")

	cpu_cores = mp.cpu_count()
	if cpu_cores > max_cpu_cores:	
		cpu_cores = max_cpu_cores;
	if cpu_cores < min_cpu_cores:	
		cpu_cores = min_cpu_cores;
	if (cpu_cores%2) and cpu_cores!=1:	
		cpu_cores -= 1
		print("[i] number cpu_cores must be even!");
	print('[cpu] %s cores available (min=%s; max=%s)' % (cpu_cores, min_cpu_cores, max_cpu_cores))

	pid = os.getpid()
	if flag_verbose > 0: 
		print("[parent#main] pid#%s " % (pid))
	print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")

	if flag_verbose not in (0,'0',False,'False','false',''):
		print('[DEBUG] level=%s' % flag_verbose)

	if prngseed in (0,'0',False,'False','false',''):
		prngseed = random.randint(1,2**32)
	random.seed(prngseed)
	if flag_verbose > -1: 
		print('[PRNGseed] %s' % prngseed)

	if flag_gmpy2:
		print('[library#] gmpy2 (full speed available)')
	else:
		print('[library#] raw python (slowly, recommend install gmpy2)')

	print('[profile#] %s' % flag_profile)

	flag_pow2bits = False
	flag_keyspace = False

	prvkey0 = False
	pubkey0 = False

	bitMin = 8
	bitMax = 120

	if len(sys.argv) > 1 :
		#bits
		try:
			pow2bits = int(sys.argv[1])
			L = 2**(pow2bits-1)
			U = 2**pow2bits
		except:
			flag_pow2bits = False
		else:
			flag_pow2bits = True
		#range
		try:
			L, U = str(sys.argv[1]).split(':')
			L = int(str(L), 16)
			U = int(str(U), 16)
			assert(len(sys.argv)>2)
			bitMin = 8
			bitMax = 256
		except:
			flag_keyspace = False
		else:
			flag_keyspace = True

		if (not flag_pow2bits) and (not flag_keyspace):
			usage()

		if U <= L:
			print("[error] 0x%x GreaterOrEqual 0x%x" % (L,U))
			usage()

		W = U - L
		try:
			Wsqrt = W**0.5
			#Wsqrt = math.sqrt(W)
			Wsqrt = int(Wsqrt)
		except:
			usage()

		# M == (L+U)/2 == L+(W/2)
		#M = (L + U)//2
		M = L + (W//2)

		if flag_pow2bits:
			pow2L = pow2bits-1
			pow2U = pow2bits
			pow2W = pow2bits-1
			print('[range] 2^%s..2^%s ; W = U - L = 0x%x (2^%s)' % (pow2L, pow2U, W, pow2W))
		if flag_keyspace:
			pow2L = int(math.log(L,2))+0
			pow2U = int(math.log(U,2))+1
			pow2W = int(math.log(W,2))+1
			pow2bits = pow2U
			print('[range] 0x%x..0x%x ; W = U - L = 0x%x (~2^%s)' % (L, U, W, pow2W))

	# without args
	else:
		flag_pow2bits = True
		flag_keyspace = False

		L = 2**(pow2bits-1)
		U = 2**pow2bits

		W = U - L
		try:
			Wsqrt = W**0.5
			#Wsqrt = math.sqrt(W)
			Wsqrt = int(Wsqrt)
		except:
			usage()

		# M == (L+U)/2 == L+(W/2)
		#M = (L + U)//2
		M = L + (W//2)

		pow2L = pow2bits-1
		pow2U = pow2bits
		pow2W = pow2bits-1
		print('[range] 2^%s..2^%s ; W = U - L = 0x%x (2^%s)' % (pow2L, pow2U, W, pow2W))


	if pow2W < bitMin or pow2W > bitMax :
		print('[error] W must be 2^%s..2^%s!' % (bitMin,bitMax))
		usage()
	if pow2W > 55 :
		print('[warn!] W = 2^%s too big! long runtime expected' % (pow2W) )

	print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")
	############
	# pre-compute set S(i) points of pow2 jumpsize

	Sp = [Gp]
	dS = [1]
	for k in xrange(255): 
		Sp.append(mul_2a(Sp[k]))
		dS.append(2*dS[k])
	print('[+] Sp-table of pow2 points - ready')

	#Sp_orig = Sp.copy(); dS_orig = dS.copy()
	Sp_orig = list(Sp); dS_orig = list(dS)
	#print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")
	############
	# timeit loop

	# delay loop checker messages of parent
	time_delay = 0.05/cpu_cores

	# delay between messages from childs
	if flag_gmpy2:
		jump_step = 10000
	else:
		jump_step = 1000

	# parent print progress
	print_eachNjumps = jump_step * cpu_cores

	starttime = time.time()
	list_sumjump, list_runtime, list_dpkgr = list(), list(), list()

	#timeit
	for n_timeit in xrange(1,Ntimeit+1):

		print("[~~~~~~~~~~~~~~~~~~~~~~[%s/%s]~~~~~~~~~~~~~~~~~~~~~]" % (n_timeit, Ntimeit))
		if flag_verbose > 1: 
			save2file('tame.txt', 'w', '')
			save2file('wild.txt', 'w', '')

		############
		# pubkey load
		if 1:

			if len(sys.argv)>2 :
				pubkey0 = str(sys.argv[2])
				print('[i] custom pubkey#%s loaded from argv2' % pow2bits)

			elif not (Ntimeit>1 and timeit_eachnewprvkey):
				try:
					pubkey0, prvkey0 = pubkeys[pow2bits]
				except:
					prvkey0 = random.randint(L,U)
					pubkey0 = getPubkey(prvkey0, True)	#   compressed
					#pubkey0 = getPubkey(prvkey0, False)	# uncompressed
					print('[i] pubkey#%s randomly generated in range [2^%s..2^%s]' % (pow2bits, pow2L, pow2U))
				else:
					print('[i] pubkey#%s loaded from default table' % pow2bits)
			else:
				if 1:
					prvkey0 = random.randint(L,U)
					pubkey0 = getPubkey(prvkey0, True)	#   compressed
					#pubkey0 = getPubkey(prvkey0, False)	# uncompressed
					print('[i] pubkey#%s randomly generated in range [2^%s..2^%s]' % (pow2bits, pow2L, pow2U))


			if prvkey0 not in (0,'0',False,'False','false',''):
				if flag_pow2bits:
					print('[prvkey#%s] 0x%064x' % (pow2bits,prvkey0))
				if flag_keyspace:
					print('[prvkey#xx] 0x%064x' % (prvkey0))

			# location in keyspace on the strip
			if flag_verbose > -1 :
				if prvkey0 not in (0,'0',False,'False','false',''):
					len100perc = 60
					size1perc = W//len100perc
					print("[i] [2^%.1f|%s%s%s|2^%.1f]" % (pow2L
						, '-'*((prvkey0-L)//size1perc)
						, 'K'
						, '-'*((U-prvkey0)//size1perc)
						, pow2U)
					);#exit(1)

			if 1:
				if flag_pow2bits:
					print('[pubkey#%s] %s' % (pow2bits,pubkey0))
				if flag_keyspace:
					print('[pubkey#xx] %s' % (pubkey0))
	
			# check format pubkey
			if len(pubkey0)==130:
				X = int(pubkey0[2:66], 16)
				Y = int(pubkey0[66:],16)
				flag_compress = False
				#print("[format] uncompressed")
			elif len(pubkey0)==66:
				X = int(pubkey0[2:66], 16)
				# calculation Y from X if pubkey is compressed
				Y = getX2Y(X,int(pubkey0[:2])-2)
				flag_compress = True
				#print("[format] compressed")
			else:
				print("[error] pubkey len(66/130) invalid!")
				usage()

			#print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")
			print("[Xcoordinate] %064x" % X)
			print("[Ycoordinate] %064x" % Y)
			#print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")

			if not is_on_curve(X,Y):
				print("[error] the given point not lies on the elliptic curve!")
				usage()
			
			# wild root
			W0p = Point(X,Y)

		############
		# profiles load

		# number kangaroos of herd T/W
		if cpu_cores == 1 or cpu_cores == 2:
			xU = xV = 1

		elif cpu_cores >= 4:
			if flag_profile == "byPollard":
				# odd int
				xU = (cpu_cores//2)-1
				xV = (cpu_cores//2)+1
				# prime int
				#xU = getPrimeInt((cpu_cores//2)-1)
				#xV = getPrimeInt((cpu_cores//2)+1)

				if flag_verbose > 0:
					print("[U] %s (0x%02x)" % (xU,xU))
					print("[V] %s (0x%02x)" % (xV,xV))

				for k in xrange(len(Sp)): 
					Sp[k] = mul_ka(xU*xV, Sp_orig[k])
					dS[k] = xU*xV*dS_orig[k]

				if flag_verbose > 0: 
					print('[+] recalc Sp-table of multiply UV')

			elif flag_profile == "NaiveSplitRange":
				xU = xV = cpu_cores//2
			elif flag_profile == "byOorschot&Wiener":
				xU = xV = cpu_cores//2
			else:
				xU = xV = cpu_cores//2


		# mean jumpsize
		if xU == xV == 1:
			# by Pollard ".. The best choice of m (mean jump size) is w^(1/2)/2 .."
			#midJsize = (Wsqrt//2)+1
			midJsize = int(round(1.0*Wsqrt/2))
		else:
			# expected of 2w^(1/2)/cores jumps
			if flag_profile == "byPollard": 
				#midJsize = int(round(1.0*((1.0*W/(xU*xV))**0.5)/2))
				midJsize = int(round(1.0*(xU+xV)*Wsqrt/4))
				#midJsize = int(round(1.0*Wsqrt/2))

			# expected of 2(w/cores)^(1/2) jumps
			elif flag_profile == "NaiveSplitRange":
				midJsize = int(round(1.0*((1.0*W/xU)**0.5)/2))

			# 
			elif flag_profile == "byOorschot&Wiener":
				midJsize = int(round(1.0*(xU+xV)*Wsqrt/4))
				#midJsize = int(round(1.0*Wsqrt/2))

			# 
			else:
				midJsize = int(round(1.0*Wsqrt/2))

		JmaxofSp = getJmaxofSp(midJsize, dS)+0

		#sizeJmax = 2**JmaxofSp
		sizeJmax = dS[JmaxofSp]

		# discriminator for filter added new distinguished points (ram economy)
		pow2dp = (pow2W//2)-2
		DPmodule = 2**pow2dp


		if flag_verbose > 0: 
			#print('[sizeJmax] 2^%s		= %s (0x%02x)' % (JmaxofSp, sizeJmax, sizeJmax))
			print('[sizeJmax] S[%s]	= %s (0x%02x)' % (JmaxofSp, sizeJmax, sizeJmax))
			print('[DPmodule] 2^%s		= %s (0x%02x)' % (pow2dp, DPmodule, DPmodule))

		############
		# create herds points

		# dT/dW - int, sum distance traveled
		dT, dW = list(), list()

		# Tp/Wp - point, sum distance traveled
		Tp, Wp = list(), list()

		# generate start points

		# Tame herd, create rand
		for k in range(xU):
			if xU == xV == 1:
				dT.append(M)
			else:
				if flag_profile == "byPollard":
					dT.append( M + (k-0)*xV )
				elif flag_profile == "NaiveSplitRange":
					dT.append( L + (W//(2*xU)) + ((k-0)*W//xU) )
				elif flag_profile == "byOorschot&Wiener":
					#dT.append( M + (1<<random.randint(1, pow2W-1))-1 )
					#dT.append( M + (1<<random.randint(1, pow2W//2))-1 )
					#dT.append( M + (k+0)*(W//(4*(xU+xV)))-1 )
					dT.append( L + (W//(2*xU)) + ((k-0)*W//xU) )
				else:
					#dT.append( (3<<(pow2bits-2)) + random.randint(1, (2**(pow2bits-1))) )	# by 57fe
					#dT.append( M + random.randint(1, W) )	# by 57fe
					dT.append( M + random.randint(1, W//2) )

			# T odd recommended (for more efficiency)
			if not (dT[k]%2):
				if flag_profile != "byPollard":	
					dT[k] += 1; 
					pass

			if flag_verbose > 1:	print('dT[%s] 0x%064x' % (k,dT[k]))

			Tp.append(mul_ka(dT[k]))
	
		# Wild herd, add rand
		for k in range(xV):
			if xU == xV == 1:
				dW.append(1)
			else:
				if flag_profile == "byPollard":
					dW.append( 1 + xU*(k-0) )
				elif flag_profile == "NaiveSplitRange":
					#dW.append( 1 + (W//(2*xU)) + ((k-0)*W//xU) )
					dW.append( 1 + ((k-0)*W//xU) )
					#dW.append( 1 )
				elif flag_profile == "byOorschot&Wiener":
					#dW.append( (1<<random.randint(1, pow2W-1))-1 )
					#dW.append( (1<<random.randint(1, pow2W//2))-1 )
					#dW.append( 1 + (k+0)*(1<<(pow2W//(xU+xV)))-1 )
					dW.append( 1 + ((k-0)*W//xU) )
				else:
					#dW.append( random.randint(1, (1<<(pow2bits-1))) )	# by 57fe
					#dW.append( random.randint(1, W) )	# by 57fe
					dW.append( random.randint(1, W//2) )

			if flag_verbose > 1:	print('dW[%s] 0x%064x' % (k,dW[k]))

			Wp.append(add_a(W0p,mul_ka(dW[k])))

		print('[+] %sT+%sW herds - ready' % (xU, xV) )

		############
		# init and clear

		# DTp/DWp - points, distinguished of Tp/Wp
		DTp, DWp = dict(), dict() # dict is hashtable of python, provides uniqueness distinguished points

		prvkey = False;
		sumjump = last_sumjump = 0 
		t0 = t1 = t2 = t1_info = t2_info = time.time()
		n_rot = 0
		repairDT = repairDW = 0
		countDT = countDW = 0

		parent_reciver = mp.Queue()
		send2childs = dict()

		############
		# run childs

		procs = list()

		## Tame herd, start childs
		for k in range(xU):
			id_uniq = k
			send2childs[id_uniq] = mp.Queue()
			proc = mp.Process(target=KANGAROO
				, name='tame-'+str(k+1)
				, args=( id_uniq 
					, Sp, dS 
					, Tp[k], dT[k]
					, DPmodule, JmaxofSp, jump_step
					, parent_reciver, send2childs[id_uniq]
				,)
			)
			proc.daemon=True
			procs.append(proc)
			proc.start()
			#proc.join()


		## Wild herd, start childs
		for k in range(xV):
			id_uniq = k+xU
			send2childs[id_uniq] = mp.Queue()
			proc = mp.Process(target=KANGAROO
				, name='wild-'+str(k+1)
				, args=( id_uniq 
					, Sp, dS 
					, Wp[k], dW[k]
					, DPmodule, JmaxofSp, jump_step
					, parent_reciver, send2childs[id_uniq]
				,)
			)
			proc.daemon=True
			procs.append(proc)
			proc.start()
			#proc.join()


		############
		# wait
		
		while(1):
			#time.sleep(0.01)
			time.sleep(time_delay)

			if not parent_reciver.empty():

				msg  = parent_reciver.get_nowait()

				id_uniq		 = msg['id_uniq']
				child_pid	 = msg['pid']
				child_name	 = msg['name']

				# send new distinguished point to parent
				#send2parent.put_nowait({'id_uniq': id_uniq, 'pid': child_pid, 'name': child_name
				#		, 'diffjumps': False, 'Xcoord': Xcoord, 'dK': dK})

				#send2parent.put_nowait(
				#		{'id_uniq': id_uniq, 'pid': child_pid, 'name': child_name
				#		, 'diffjumps': countj-last_countj, 'Xcoord': False, 'dK': False}

				if msg['diffjumps'] != False:
					sumjump	+= msg['diffjumps']
					if flag_verbose > 1: 
						print("\n[childs][%s#%s] +%s jumps %40s" % (child_name,child_pid, msg['diffjumps'], ' '));
				else:
					Xcoord	 = msg['Xcoord']
					dK	 = msg['dK']

					if flag_verbose > 1: 
						print("\n[childs][%s#%s] newDP: X=%064x " % (child_name,child_pid, Xcoord));
						print("\n[childs][%s#%s] newDP: dK=0x%x " % (child_name,child_pid, dK));

					# Tame herd
					if id_uniq < xU:
						try:
							DTp[Xcoord]
						except:
							# add new distinguished point
							DTp[Xcoord] = dK
							if flag_verbose > 1: 
								save2file('tame.txt', 'a', '%064x %s\n'%(Xcoord,dK) )
						else:
							# repeat detected!
							repairDT += 1
							if flag_verbose > 0: 
								#print('\n[parent][%s#%s] repair#%s/%s: 0x%064x ' % (child_name,child_pid, repairDT,repairDW, Xcoord));
								print('\n[parent][%s#%s] repair#%s/%s: 0x%x ' % (child_name,child_pid, repairDT,repairDW, dK));
							# need fix same sequences
							send2childs[id_uniq].put_nowait({'dK': dK})

					# Wild herd
					#if id_uniq >= xU:
					else:
						try:
							DWp[Xcoord]
						except:
							# add new distinguished point
							DWp[Xcoord] = dK
							if flag_verbose > 1: 
								save2file('wild.txt', 'a', '%064x %s\n'%(Xcoord,dK) )
						else:
							# repeat detected!
							repairDW += 1
							if flag_verbose > 0: 
								#print('\n[parent][%s#%s] repair#%s/%s: 0x%064x ' % (child_name,child_pid, repairDT,repairDW, Xcoord));
								print('\n[parent][%s#%s] repair#%s/%s: 0x%x ' % (child_name,child_pid, repairDT,repairDW, dK));
							# need fix same sequences
							send2childs[id_uniq].put_nowait({'dK': dK})


					# compare distinguished points, Tame herd & Wild herd
					compare = list(set(DTp) & set(DWp))
					if len(compare) > 0: 
						dDT = DTp[compare[0]]
						dDW = DWp[compare[0]]
						if	dDT > dDW:
							prvkey = dDT - dDW
						elif	dDW > dDT:
							prvkey = dDW - dDT
						else:
							print("\n[FATAL_ERROR] dDW == dDT !!! (0x%x)\n" % (dDW));exit(-1)

			#exit(1)

			# end
			#if prvkey:
			#	runtime = time.time()-t0 
			#	break

			#exit(1)
			
			#if not (sumjump % print_eachNjumps) and sumjump != 0:
			if (sumjump-last_sumjump)>=print_eachNjumps  and  sumjump != 0:

				t2 = t2_info = time.time()

				countDT = len(DTp); countDW = len(DWp);

				# info
				if (flag_verbose > 0 and (t2_info-t1_info)>10.0)  or prvkey:
					printstr  = '\r[i] DP %sT+%sW=%s+%s=%s; dp/kgr=%.1f' % ( xU, xV
							, countDT,countDW, countDT+countDW, 1.0*(countDT+countDW)/2
							)
					printstr += ' '*50
					print(printstr)
					t1_info = t2_info

				# indicator, progress, time
				if (t2-t1)>=1.0  or prvkey:
					printstr  = '\r'

					if	n_rot==1:	printstr += "[\]"; n_rot=2;
					elif	n_rot==2:	printstr += "[|]"; n_rot=3;
					elif	n_rot==3:	printstr += "[/]"; n_rot=0;
					else:			printstr += "[-]"; n_rot=1;

					printstr += '[%s ' % ( time_format(t2-t0, (1,1,1,1,1,1,0,0)) )
					if t2-t1 != 0:
						printstr += '; %s j/s' % prefSI((sumjump-last_sumjump)/(t2-t1))
					else:
						printstr += '; %s j/s' % prefSI((sumjump-last_sumjump)/1)

					#printstr += '; %sj of %sj %.1f%%' % (
					printstr += '; %sj %.1f%%' % (
							 sumjump if sumjump<10**3 else prefSI(sumjump)
							#, 2*Wsqrt if 2*Wsqrt < 10**3 else prefSI(2*Wsqrt)
							, (1.0*sumjump/(2*Wsqrt))*100
							)
					if 1 or flag_verbose < 1: 
						printstr += '; dp/kgr=%.1f' % ( 1.0*(countDT+countDW)/2 )
					#printstr += 'lost_TIME_left'
					timeleft = (t2-t0)*(1-(1.0*sumjump/(2*Wsqrt)))/(1.0*sumjump/(2*Wsqrt))
					if timeleft > 0:
						printstr += ';%s ]  ' % ( time_format(timeleft, (1,1,1,1,1,1,0,0)) )
					else:
						printstr += ';%s ]  ' % ( time_format(0, (1,1,1,1,1,1,0,0)) )
					if sys.version_info[0] == 2:
						print(printstr, end='')
						sys.stdout.flush()
					else:
						print(printstr, end=''
						, flush=True )
					t1 = t2
					last_sumjump = sumjump


			# end
			if prvkey:
				runtime = time.time()-t0 
				break


		# kill childs
		for proc in procs:
			if proc.is_alive(): proc.terminate()



		if sumjump == 0:	sumjump = print_eachNjumps
		countDT = len(DTp); countDW = len(DWp);

		# save stat for avg
		list_sumjump.append(sumjump)
		list_runtime.append(runtime)
		list_dpkgr.append(1.0*(countDT+countDW)/2)
		
		print('')
		print('[prvkey#%s] 0x%064x' % (pow2bits,prvkey) )

		# location in keyspace on the strip
		if 1:
			if flag_verbose > -1:
				len100perc = 60
				size1perc = W//len100perc
				print("[i] [2^%.1f|%s%s%s|2^%.1f]" % (
					 math.log(L,2)
					, '-'*((prvkey-L)//size1perc)
					, 'K'
					, '-'*((U-prvkey)//size1perc)
					, math.log(U,2)
					)
				);#exit(1)
		

		# double-check privkey
		if prvkey0 not in (0,'0',False,'False','false',''):
			if prvkey != prvkey0:
				print('[origin#%s] 0x%064x' % (pow2bits,prvkey0))
				print('[prvkey-check] failed!')

		# double-check pubkey
		pubkey = getPubkey(prvkey,flag_compress)
		if pubkey != pubkey0:
			print('[pubkey#%s] %s' % (pow2bits,pubkey))
			print('[origin#%s] %s' % (pow2bits,pubkey0))
			print('[pubkey-check] failed!')
		#else:
		if 1:
			save2file('results.txt', 'a', ('%064x:%s\n' % (prvkey, pubkey), '---------------\n'))

		# finish stat
		printstr = '[i] %s j/s; %sj of %sj %.1f%%; DP T+W=%s+%s=%s; dp/kgr=%.1f' % (
			 prefSI(sumjump/1) if runtime==0 else prefSI(sumjump/runtime)
			, sumjump if sumjump<10**3 else prefSI(sumjump)
			, 2*Wsqrt if 2*Wsqrt < 10**3 else prefSI(2*Wsqrt)
			, (1.0*sumjump/(2*Wsqrt))*100
			, countDT,countDW, countDT+countDW, 1.0*(countDT+countDW)/2
			)
		printstr += '  '
		print(printstr)
		#print('[runtime]%s' % time_format(runtime))
		print('[runtime]%s' % time_format(runtime, (1,1,1,1,1,1,2,2)))


	print("[################################################]")

	#avgTime = (time.time()-starttime)/Ntimeit
	avgTime = 1.0 * sum(runtime for runtime in list_runtime) / len(list_runtime)
	avgJump = 1.0 * sum(sumjump for sumjump in list_sumjump) / len(list_sumjump)
	avgDPkg = 1.0 * sum(rundpkg for rundpkg in list_dpkgr) / len(list_dpkgr)
	#D = sum((xi - avgJump) ** 2 for xi in list_sumjump)*1.0 / len(list_sumjump)

	if Ntimeit > 1:
		##print('[(avg)jump] %.0f' % (avgJump) )
		#print('[(avg)jump] %s ' % (int(avgJump) if avgJump<10**3 else prefSI(avgJump)) )
		##print('[(avg)jum2] %.1f +/- %.1f' % (avgJump, (D/(len(list_sumjump)-1))**0.5) )
		#print('[(avg)dpkg] %s ' % (int(avgDPkg) if avgDPkg<10**3 else prefSI(avgDPkg)) )
		#print('[(avg)time]%s' % time_format(avgTime, (1,1,1,1,1,1,2,2)) )

		print("[averages] expected of 2w^(1/2) group operations")
		print("-------|--------/--------|---------------------------------/---------------------------------|")
		print("   W   |jump avg/2w^(1/2)| time                         avg/2w^(1/2)                         |")
		print("-------|--------/--------|---------------------------------/---------------------------------|")
		if 1:
			i = pow2W
			xi = 1
			print('%s2^%03d |  %06s/ %06s |%032s /%032s |' % 
					(	'>' if i==pow2W else ' '
						,i
						,int(avgJump) if int(avgJump*xi)<10**3 else prefSI(avgJump*xi)
						,int(2*(2**i)**0.5) if int(2*(2**i)**0.5)<10**3 else prefSI(2*(2**i)**0.5)
						,time_format( avgTime * xi , (1,1,1,1,1,1,2,2)) 
						,time_format( (avgTime * xi)/(avgJump/(2*Wsqrt)) , (1,1,1,1,1,1,2,2)) 
					) 
			)
		print("----------------------------------------------------------------------------------------------")
	else:
		pass

	print("[################################################]")
	print('[date] %s'%(time.ctime()))
	print("[################################################]")

	if 1:
		#bitMin = 10
		#bitMax = 256

		try:
			print('');raw_input('Press ENTER to get [prognose] or Ctrl+C to [exit] ...');print('')
		except:
			print('\n[exit] exit')
			exit(0)

		print("[prognose] expected of 2w^(1/2) group operations")

		print("-------|--------/--------|---------------------------------/---------------------------------|")
		print("   W   |jump avg/2w^(1/2)| time                         avg/2w^(1/2)                         |")
		print("-------|--------/--------|---------------------------------/---------------------------------|")
		for i in xrange(bitMin,pow2W):
			xi = ((2**i)**0.5) / Wsqrt
			print('%s2^%03d |  %06s/ %06s |%032s /%032s |' % 
					(	' '
						,i
						,int(avgJump*xi) if int(avgJump*xi)<10**3 else prefSI(avgJump*xi)
						,int(2*(2**i)**0.5) if int(2*(2**i)**0.5)<10**3 else prefSI(2*(2**i)**0.5)
						,time_format( avgTime * xi , (1,1,1,1,1,1,2,2)) 
						,time_format( (avgTime * xi)/(avgJump/(2*Wsqrt)) , (1,1,1,1,1,1,2,2)) 
					) 
			)
		for i in xrange(pow2W,bitMax+1):
			xi = ((2**i)**0.5) / Wsqrt
			print('%s2^%03d |  %06s/ %06s |%032s /%032s |' % 
					(	'>' if i==pow2W else ' '
						,i
						,int(avgJump*xi) if int(avgJump*xi)<10**3 else prefSI(avgJump*xi)
						,int(2*(2**i)**0.5) if int(2*(2**i)**0.5)<10**3 else prefSI(2*(2**i)**0.5)
						,time_format( avgTime * xi , (1,1,1,1,1,1,2,2)) 
						,time_format( (avgTime * xi)/(avgJump/(2*Wsqrt)) , (1,1,1,1,1,1,2,2)) 
					) 
			)
		print("----------------------------------------------------------------------------------------------")

	print("[################################################]")
	#print('[date] %s'%(time.ctime()))
	print('[exit] exit')
	#print('');raw_input('Press ENTER to continue...');print('')
	exit(0)
