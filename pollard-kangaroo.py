#!/usr/bin/python

# based on code by 57fe, 2019
# fe57.org/forum/thread.php?board=4&thema=1#1

#######################
# print() compatibility python 2/3
from __future__ import print_function
#######################
# users settings

pow2bits	= 32	# bits/suborder/exp key

flag_profile	= "optimal"	# settings by Telariust ; x0.9 expected of 2w^(1/2) group operations
#flag_profile	= "by 57fe"	# settings by 57fe ; x1.3 expected of 2w^(1/2) group operations
#flag_profile	= "standart"	# settings by Pollard, Oorschot, Wiener ; x3.0 expected of 2w^(1/2) group operations


#######################
# service setings

Ntimeit		= 10		# times for avg runtime
timeit_eachnewprvkey = True	# gen new privkey each loop?

prngseed	= 0	# 0 for random, or any for replay results
flag_debug	= 0	# 0, 1, 2

version = '0.82'

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

# coincurve is good
# [windows>python.exe - m] pip install coincurve
try:
	from coincurve import PrivateKey, PublicKey
	from coincurve.utils import int_to_bytes, hex_to_bytes, bytes_to_int, bytes_to_hex, int_to_bytes_padded
except:
	flag_coincurve = False
	if (not flag_gmpy2):
		print("[warn] lib coincurve not found. full speed is not achievable!")
else:
	flag_coincurve = True

# boosted coincurve
try:
	from cffi import FFI
	ffi = FFI()
except:
	flag_cffi = False
	if (not flag_gmpy2) and flag_coincurve:
		print("[warn] lib cffi not found. full speed is not achievable!")
else:
	flag_cffi = True

# debug lib
if 0:
	flag_cffi = 0
	flag_coincurve = 0
	flag_gmpy2 = 0

if 0:
	from multiprocessing import Pool
	from multiprocessing import cpu_count
	from multiprocessing import freeze_support


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
# secp256k1

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

if (not flag_gmpy2) and flag_coincurve:
	Gp = PublicKey.from_point(Gx, Gy) # basePoint 
	#Gp = PublicKey.from_valid_secret(int_to_bytes_padded(1)) # basePoint 
else:
	Gp = Point(Gx,Gy) 
	Zp = Point(0,0)	# zero-point, infinite in real x,y - plane


#######################
# math, raw python

# return (g, x, y) a*x + b*y = gcd(x, y)
def egcd(a, b):
	if a == 0:
		return (b, 0, 1)
	else:
		g, x, y = egcd(b % a, a)
		return (g, y - (b // a) * x, x)

# from 57fe origin script
def invert_old(b, p=modulo):
	while b < 0:
		b += p
	g, x, _ = egcd(b, p)
	if g == 1:
		return x % p


# from arulberoECC library
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
# ..so option "0S" (without pow2) more preferred


# A + A -> A (1I, 2M, 1S)
# A + A -> A (1I, 3M, 0S)
def add_a(A, B, p=modulo):
	R = Point()
	dx = B.x - A.x
	dy = B.y - A.y	
	if flag_gmpy2:
		c = dy * gmpy2.invert(dx, p) % p	# 1I,1M
	else:
		#c = dy * invert_old(dx, p) % p		# old invert
		c = dy * invert(dx, p) % p		# 1I,1M

	#R.x = (c**2 - A.x - B.x) % p	# 0M,1S
	R.x = (c*c - A.x - B.x) % p	# 1M,0S
	#R.x = (int(math.pow(c,2)) - A.x - B.x) % p	# slow

	R.y = (c*(A.x - R.x) - A.y) % p	# 1M
	return R


# 2 * A -> A (1I, 5M, 2S)
# 2 * A -> A (1I, 7M, 0S)
def mul_2a(A, p=modulo):
	R = Point()
	if flag_gmpy2:
		#c = 3 * A.x**2 * gmpy2.invert(2*A.y, p) % p	# 1I,3M,1S
		c = 3 * A.x * A.x * gmpy2.invert(2*A.y, p) % p	# 1I,4M,0S
	else:
		#c = 3 * A.x**2 * invert_old(2*A.y, p) % p	# old invert
		#c = 3 * A.x**2 * invert(2*A.y, p) % p		# 1I,3M,1S
		c = 3 * A.x * A.x * invert(2*A.y, p) % p	# 1I,4M,0S

	#R.x = (c**2 - 2*A.x) % p	# 1M,1S
	R.x = (c*c - 2*A.x) % p		# 2M,0S

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
	mcsec= int((((time%1)*1000)%1)*1000)
	res  = ''	
	if v[0]: 
		y_tmp = (sec//(60*60*24*30))//12
		if y_tmp>0: res += ' '+'%s'%(y_tmp if y_tmp<10**3 else prefSI(y_tmp))	+'y'	# year
	if v[1]:
		m_tmp = (sec//(60*60*24*30))%12
		if m_tmp>0 or y_tmp>0: res += ' '+'%02s'%str(m_tmp)			+'m'	# month
	if v[2]: res += ' '+'%02s'%str((sec//(60*60*24))%30)				+'d'	# day
	if v[3]: res += ' '+'%02d'%int((sec//(60*60))%24)				+''	# hour
	if v[4]: res += ':'+'%02d'%int((sec//(60*1))%60)				+''	# min
	if v[5]: res += ':'+'%02d'%int((sec//(1*1))%60)					+'s'	# sec
	if v[6]: res += ' '+'%03d'%msec							+'ms'	# msec
	if v[7]: res += ' '+'%03d'%mcsec						+'mcs'	# mcsec
	return res
#print('[time] %s'%time_format(int(sys.argv[1])));exit(1)


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
	elif flag_coincurve:
		if flag_cffi:
			tmp_pubkey = ffi.buffer(itpoint.public_key, 64)[:]
			Xcoord = bytes_to_int(tmp_pubkey[31::-1])
			#Ycoord = bytes_to_int(tmp_pubkey[:31:-1])
		else:
			Xcoord, Ycoord = itpoint.point()
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
	elif flag_coincurve:
		Ptmp = PublicKey.from_valid_secret(int_to_bytes_padded(new_prvkey))
		if flag_cffi:
			tmp_pubkey = ffi.buffer(Ptmp.public_key, 64)[:]
			Xcoord = bytes_to_int(tmp_pubkey[31::-1]); 
			Ycoord = bytes_to_int(tmp_pubkey[:31:-1]);
		else:
			Xcoord, Ycoord = Ptmp.point()
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




#######################
# KANGAROOS

def KANGAROOS():

	# settings by Telariust
	# x0.9 expected of 2w^(1/2) group operations
	if flag_profile == "optimal":
		# pow2 size herd T+W (number of kangaroos in T/W herd), affects max size jump, affects discriminator
		pow2kang = 0

		# HTmax,HWmax - max number of kangaroos in T/W herd
		HTmax = HWmax = 2**pow2kang

		# mean jumpsize
		# from Pollard ".. The best choice of m (jump size) is (w^(1/2))/2 .."
		#midJsize = Wsqrt/2 # m = w^(1/2)/2
		# for N kangaroos
		#midJsize = (HTmax+HWmax)*Wsqrt/4	# with N cpu is m = N(w^(1/2))/4

		# max jumpsize is 2m
		#maxJsize = int(2*midJsize)
		maxJsize = Wsqrt*8	# x8 the fastest! empirically

		# pow max jumpsize
		pow2Jmax = int(math.log(maxJsize,2))+1

		range1 = M	# standart
		range2 = 1	# 1 for T1+W1

		# discriminator for filter added new distinguished points (ram economy)
		pow2dp = ((pow2W - 2*pow2kang)//2)-2	# 
		DP_rarity = 2**pow2dp

	# settings by 57fe
	# x1.3 expected of 2w^(1/2) group operations
	elif flag_profile ==  "by 57fe":
		# pow2 size herd T+W (number of kangaroos in T/W herd), affects max size jump, affects discriminator
		pow2kang = 3

		# HTmax,HWmax - max number of kangaroos in T/W herd
		HTmax = HWmax = 2**pow2kang

		# powmax jumpsize
		pow2Jmax = (pow2W//2) + pow2kang	# by 57fe
		# max jumpsize is 2m
		maxJsize = 2**pow2Jmax

		range1 = M	# standart
		range2 = W	# by 57fe

		# discriminator for filter added new distinguished points (ram economy)
		pow2dp = ((pow2W - 2*pow2kang)//2)-2	# by 57fe
		DP_rarity = 2**pow2dp

	# settings by Pollard, Oorschot, Wiener
	# x3 expected of 2w^(1/2) group operations
	#elif flag_profile ==  "standart":
	else:
		# pow2 size herd T+W (number of kangaroos in T/W herd), affects max size jump, affects discriminator
		pow2kang = 0

		# HTmax,HWmax - max number of kangaroos in T/W herd
		HTmax = HWmax = 2**pow2kang

		# mean jumpsize
		# from Pollard ".. The best choice of m (jump size) is (w^(1/2))/2 .."
		#midJsize = Wsqrt/2 # m = w^(1/2)/2
		# for N kangaroos
		midJsize = (HTmax+HWmax)*Wsqrt/4	# with N cpu is m = N(w^(1/2))/4

		# max jumpsize is 2m
		maxJsize = int(2*midJsize)

		# pow max jumpsize
		pow2Jmax = int(math.log(maxJsize,2))+1

		range1 = M	# standart
		range2 = 1	# 1 for T1+W1

		# discriminator for filter added new distinguished points (ram economy)
		pow2dp = ((pow2W - 2*pow2kang)//2)-2	# 
		DP_rarity = 2**pow2dp


	if flag_debug > 0: 
		print('[kangaroos] 2^%s = %s in Tame/Wild herd' % (pow2kang, 2**pow2kang))
		if flag_debug > 0: 
			print('[DP_rarity] 2^%s = %s' % (pow2dp, DP_rarity))
			print('[max_Jsize] 2^%s = %s' % (pow2Jmax, maxJsize))

	# dt/dw - int, current jumpsize
	# dT/dW - int, sum dt/dw as int, distance traveled
	# Tp/Wp - point, sum dt/dw as points, distance traveled
	Tp, dT, dt = list(), list(), list()
	Wp, dW, dw = list(), list(), list()

	# generate random start points
	if flag_pow2bits:
		#if flag_debug > 1:	print('dT[k] (3/4)*(2^bits) + rng((1/2)*(2^bits))')	# by 57fe
		#if flag_debug > 1:	print('dT[k] (3/4)*(2^bits) + rng(?)')	# M + rng(?)
		pass
	if flag_keyspace:
		# M == L+(W/2) == (L+U)/2
		#if flag_debug > 1:	print('dT[k] M + rng(1, W)')	# by 57fe
		#if flag_debug > 1:	print('dT[k] M + rng(?)')	# M + rng(?)
		pass

	# Tame herd
	for k in xrange(HTmax):
		if flag_pow2bits:
			#dT.append((3<<(pow2bits-2)) + random.randint(1,(2**(pow2bits-1))))	# by 57fe
			dT.append( range1 + random.randint(1, range2))	# ? + ?
		if flag_keyspace:
			#dT.append( M + random.randint(1, W))	# by 57fe
			dT.append( range1 + random.randint(1, range2))	# ? + ?
		if (not flag_gmpy2) and flag_coincurve:
			Tp.append(Gp.multiply(int_to_bytes(dT[k])))
		else:
			Tp.append(mul_ka(dT[k]))

		dt.append(0)

		#if flag_debug > 1:	print('dT[%s] 0x%x + rng(1,0x%x) = 0x%x' % (k+1, 3<<(pow2bits-2), 2**(pow2bits-1), dT[k]))	# by 57fe
		#if flag_debug > 1:	print('dT[%s] 0x%x + rng(1,0x%x) = 0x%x' % (k+1, M, W, dT[k]))	# by 57fe
		if flag_debug > 1:	print('dT[%s] 0x%064x' % (k+1, dT[k]))

		#location in keyspace on the strip
		if flag_debug > 0 and (not (HTmax==1 and HWmax==1)):	
				len100perc = 56
				size1perc = W//len100perc
				percKey = (dT[k]-M)//size1perc
				print('dT[%s] [(3/4)2^%s|%s%s%s|2^%s|%s%s%s|(5/4)2^%s]' % (k+1
					, pow2U
					, '-'*( percKey if percKey in range(0,len100perc//2) else len100perc//2 )
					, ( 'T' if percKey in range(0,len100perc//2) else '-' )
					, '-'*( (len100perc//2)-percKey if percKey in range(0,len100perc//2) else 0 )
					, pow2U
					, '-'*( percKey-(len100perc//2) if percKey in range(len100perc//2,len100perc) else len100perc//2 )
					, ( 'T' if percKey in range(len100perc//2,len100perc) else '-' )
					, '-'*( len100perc-percKey if percKey in range(len100perc//2,len100perc) else 0 )
					, pow2U
					)
				);#exit(1)
	# Wild herd
	for k in xrange(HWmax):
		if flag_pow2bits:
			#dW.append(random.randint(1, (1<<(pow2bits-1))))	# by 57fe
			dW.append(random.randint(1, range2))	# ? + ?
		if flag_keyspace:
			#dW.append(random.randint(1, W))			# by 57fe
			dW.append(random.randint(1, range2))	# ? + ?

		if (not flag_gmpy2) and flag_coincurve:
			Wp.append(PublicKey.combine_keys([W0p,Gp.multiply(int_to_bytes(dW[k]))]))
		else:
			Wp.append(add_a(W0p,mul_ka(dW[k])))

		dw.append(0)

		#if flag_debug > 1:	print('dW[%s] 0x%x + rng(1,0x%x) = 0x%x' % (k+1, 3<<(pow2bits-2), 2**(pow2bits-1), dW[k]))	# by 57fe
		#if flag_debug > 1:	print('dW[%s] 0x%x + rng(1,0x%x) = 0x%x' % (k+1, M, W, dW[k]))	# by 57fe
		if flag_debug > 1:	print('dW[%s] 0x%064x' % (k+1, dW[k]))


	#print('[kangaroos] 2^%s = %s in Tame/Wild herd' % (pow2kang, 2**pow2kang))
	print('[+] T%s+W%s herds - ready' % (HTmax,HWmax))
	#exit(1)

	# DTp/DWp - points, distinguished of Tp/Wp
	DTp, DWp = dict(), dict() # dict is hashtable of python, provides uniqueness distinguished points

	t0 = t1 = t2 = t1_info = t2_info = time.time()
	n_jump = last_jump = n_loop = 0
	prvkey = False;

	pow2repair = 0 # for fix same sequences
	nTrepair = nWrepair = 0

	# main loop
	while (1):
		n_loop += 1
		if flag_debug > 2: 
			print('\r[debug] %s loops; %s jumps' % (n_loop, n_jump))
		
		# Tame herd
		for k in xrange(HTmax):
			if flag_debug > 2: print('\r[debug] T%s=%s, %s repairs' % (HTmax,k+1,nTrepair))
			n_jump += 1

			# Xcoord
			Xcoord = getXcoord(Tp[k])

			pw = Xcoord % pow2Jmax
			pw = int(pw)
			dt[k] = 1<<pw

			# check, is it distinguished point?
			if Xcoord % DP_rarity == 0:
				# uniqueness?
				while(1):
					try:
						DTp[Xcoord]
					except:
						break
					else:
						# repeat detected!
						nTrepair += 1
						if flag_debug > 0: 
							print('\r[tame#%s] repair: 0x%064x' % (k+1,Xcoord));
						# need fix same sequences
						dT[k] += 1<<pow2repair
						
						if (not flag_gmpy2) and flag_coincurve:
							Tp[k] = PublicKey.combine_keys([Sp[pow2repair], Tp[k]])
						else:
							Tp[k] = add_a(Sp[pow2repair], Tp[k])

						# Xcoord
						Xcoord = getXcoord(Tp[k])

						pw = Xcoord % pow2Jmax
						pw = int(pw)
						dt[k] = 1<<pw

				# add new distinguished point
				DTp[Xcoord] = dT[k]

				if flag_debug > 1: 
					printstr  = '\r[tame] T%s/W%s=%s/%s' % (HTmax,HWmax, len(DTp),len(DWp))
					printstr += '' if (HTmax==1 and HWmax==1) else '; %s/%s repairs' % (nTrepair,nWrepair)
					printstr += '; %064x 0x%x' % (Xcoord,dT[k])
					print(printstr)
					save2file('tame.txt', 'a', '%064x %s\n'%(Xcoord,dT[k]) )
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
						print("\r[error] dDW == dDT !!! (0x%x)"%dDW);exit(-1)

			if prvkey: break
			dT[k] += dt[k]
			if (not flag_gmpy2) and flag_coincurve:
				Tp[k] = PublicKey.combine_keys([Sp[pw], Tp[k]])
			else:
				Tp[k] = add_a(Sp[pw], Tp[k])
		if prvkey: break
			
		# Wild herd
		for k in xrange(HWmax):
			if flag_debug > 2: print('\r[debug] W%s=%s, %s repairs'%(HWmax,k+1,nWrepair))
			n_jump += 1

			# Xcoord
			Xcoord = getXcoord(Wp[k])

			pw = Xcoord % pow2Jmax
			pw = int(pw)
			dw[k] = 1<<pw

			# add new distinguished point
			if Xcoord % DP_rarity == 0:
				# uniqueness?
				while(1):
					try:
						DWp[Xcoord]
					except:
						break
					else:
						# repeat detected!
						nWrepair += 1
						if flag_debug > 0: 
							print('\r[wild#%s] repair: 0x%064x' % (k+1,Xcoord));
						# need fix same sequences
						dW[k] += 1<<pow2repair

						if (not flag_gmpy2) and flag_coincurve:
							Wp[k] = PublicKey.combine_keys([Sp[pow2repair], Wp[k]])
						else:
							Wp[k] = add_a(Sp[pow2repair], Wp[k])
						
						# Xcoord
						Xcoord = getXcoord(Wp[k])

						pw = Xcoord % pow2Jmax
						pw = int(pw)
						dw[k] = 1<<pw

				# add new distinguished point
				DWp[Xcoord] = dW[k]

				if flag_debug > 1: 
					printstr  = '\r[wild] T%s/W%s=%s/%s' % (HTmax,HWmax, len(DTp),len(DWp))
					printstr += '' if (HTmax==1 and HWmax==1) else '; %s/%s repairs' % (nTrepair,nWrepair)
					printstr += '; %064x 0x%x' % (Xcoord,dW[k])
					print(printstr)
					save2file('wild.txt', 'a', '%064x %s\n'%(Xcoord,dW[k]) )
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
						print("\r[error] dDW == dDT !!! (0x%x)"%dDW);exit(-1)

			if prvkey: break
			dW[k] += dw[k]
			if (not flag_gmpy2) and flag_coincurve:
				Wp[k] = PublicKey.combine_keys([Sp[pw], Wp[k]])
			else:
				Wp[k] = add_a(Sp[pw], Wp[k])
		#if prvkey: break

		# info
		t2 = t2_info = time.time()
		if (flag_debug > 0 and (t2_info-t1_info)>10)  or prvkey:
			printstr  = '\r[i] DP T%s+W%s=%s+%s=%s; dp/kgr=%.1f' % (
					 HTmax,HWmax, len(DTp),len(DWp), len(DTp)+len(DWp), (len(DTp)+len(DWp))/(HTmax+HWmax)
					)
			printstr += ' '*60 if (HTmax==1 and HWmax==1) else '; %s/%s repairs %s' % (nTrepair,nWrepair,' '*45)
			print(printstr)
			t1_info = t2_info

		# indicator, progress, time
		#t2 = time.time()
		if (t2-t1)>1  or prvkey:
			printstr  = '\r[~] %s j/s' % prefSI((n_jump-last_jump)/(t2-t1))
			#printstr += '; %sj of %sj %.1f%%' % (
			printstr += '; %sj %.1f%%' % (
					 n_jump if n_jump<10**3 else prefSI(n_jump)
					#, 2*Wsqrt if 2*Wsqrt < 10**3 else prefSI(2*Wsqrt)
					, (1.0*n_jump/(2*Wsqrt))*100
					)
			if 1 or flag_debug < 1: 
				printstr += '; dp/kgr=%.1f' % ( 1.0*(len(DTp)+len(DWp))/(HTmax+HWmax) )
			printstr += '; [%s ' % ( time_format(t2-t0, (1,1,1,1,1,1,0,0)) )
			printstr += 'lost_TIME_left'
			timeleft = (t2-t0)*(1-(1.0*n_jump/(2*Wsqrt)))/(1.0*n_jump/(2*Wsqrt))
			if timeleft > 0:
				printstr += '%s ]  ' % ( time_format(timeleft, (1,1,1,1,1,1,0,0)) )
			else:
				printstr += '%s ]  ' % ( time_format(0, (1,1,1,1,1,1,0,0)) )
			if sys.version_info[0] == 2:
				print(printstr, end='')
				sys.stdout.flush()
			else:
				print(printstr, end=''
				, flush=True )
			t1 = t2
			last_jump = n_jump

		if prvkey: break

	return prvkey, n_jump, (time.time()-t0), len(DTp),len(DWp), HTmax,HWmax, nTrepair,nWrepair


#######################
#main

if __name__ == '__main__':

	#print('[os] %s' % os.name)
	if os.name == 'nt':
		#freeze_support()
		pass
                                       ##
	print("[################################################]")
	print("[#    Pollard-kangaroo PrivKey Recovery Tool    #]")
	print("[#            bitcoin ecdsa secp256k1           #]")
	#print("[#          based on code by 57fe 2019          #]")
	print("[#                  singlecore                  #]");
	#print("[#                  multicore                   #]");
	print("[#                    ver%04s                   #]"%version);
	print("[################################################]")

	if len(sys.argv) > 1 and str(sys.argv[1]) in ('--help','-h','/?') :
		usage()

	print('[date] {}'.format(time.ctime()))
	print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")

	if flag_debug not in (0,'0',False,'False','false',''):
		print('[DEBUG] level=%s' % flag_debug)

	if prngseed in (0,'0',False,'False','false',''):
		prngseed = random.randint(1,2**32)
	random.seed(prngseed)
	if flag_debug > -1: 
		print('[PRNGseed] %s' % prngseed)

	if flag_gmpy2:
		print('[library#] gmpy2 (full speed available)')
	elif flag_coincurve and (not flag_cffi):
		print('[library#] coincurve without cffi (fast, but gmpy2 faster)')
	elif flag_coincurve and flag_cffi:
		print('[library#] coincurve with cffi (fast, but gmpy2 faster)')
	else:
		print('[library#] raw python (slowly, recommend install gmpy2 or coincurve)')

	print('[profile#] %s'%flag_profile)

	flag_pow2bits = False
	flag_keyspace = False

	prvkey0 = False
	pubkey0 = False

	bitMin = 10
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
			bitMin = 10
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
	starttime = time.time()

	Sp = [Gp]
	for k in xrange(255): 
		if (not flag_gmpy2) and flag_coincurve:
			Sp.append(Sp[k].multiply(int_to_bytes(2)))
		else:
			Sp.append(mul_2a(Sp[k]))
	print('[+] Sp-table of pow2 points - ready')

	#print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")

	list_runjump, list_runtime, list_dpkgr = list(), list(), list()

	#timeit
	for i in xrange(Ntimeit):

		print("[~~~~~~~~~~~~~~~~~~~~~~[%s/%s]~~~~~~~~~~~~~~~~~~~~~]"%(i+1,Ntimeit))
		if flag_debug > 1: 
			save2file('tame.txt', 'w', '')
			save2file('wild.txt', 'w', '')

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

			# location in keyspace on the strip
			if flag_debug > -1 :
				if prvkey0 not in (0,'0',False,'False','false',''):
					len100perc = 60
					size1perc = W//len100perc
					print("[i] [2^%.1f|%s%s%s|2^%.1f]" % (pow2L
						, '-'*((prvkey0-L)//size1perc)
						, 'K'
						, '-'*((U-prvkey0)//size1perc)
						, pow2U)
					);#exit(1)

			if flag_pow2bits:
				if prvkey0 not in (0,'0',False,'False','false',''):
					print('[prvkey#%s] 0x%064x' % (pow2bits,prvkey0))
				print('[pubkey#%s] %s' % (pow2bits,pubkey0))
			if flag_keyspace:
				if prvkey0 not in (0,'0',False,'False','false',''):
					print('[prvkey#xx] 0x%064x' % (prvkey0))
				print('[pubkey#xx] %s' % (pubkey0))
	
			# check format pubkey
			if len(pubkey0)==130:
				X = int(pubkey0[2:66], 16)
				Y = int(pubkey0[66:],16)
				flag_compress = False
				print("[format] uncompressed")
			elif len(pubkey0)==66:
				X = int(pubkey0[2:66], 16)
				# calculation Y from X if pubkey is compressed
				Y = getX2Y(X,int(pubkey0[:2])-2)
				flag_compress = True
				print("[format] compressed")
			else:
				print("[error] pubkey len(66/130) invalid!")
				usage()

			#print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")
			print("[Xcoordinate] %064x" % X)
			print("[Ycoordinate] %064x" % Y)
			#print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")

			# wild root
			if (not flag_gmpy2) and flag_coincurve:
				W0p = PublicKey.from_point(X,Y)
			else:
				W0p = Point(X,Y)


		# call KANGAROOS()
		prvkey, runjump, runtime, lenT,lenW, HTmax,HWmax, nTrepair,nWrepair = KANGAROOS()

		# save stat for avg
		list_runjump.append(runjump)
		list_runtime.append(runtime)
		list_dpkgr.append(1.0*(lenT+lenW)/(HTmax+HWmax))
		
		print('')
		print('[prvkey#%s] 0x%064x' % (pow2bits,prvkey) )
		save2file('results.txt', 'a', ('%064x\n'%prvkey, '---------------\n'))

		# double-check privkey
		if prvkey0 not in (0,'0',False,'False','false',''):
			if prvkey != prvkey0:
				print('[origin#%s] 0x%064x' % (pow2bits,prvkey0))
				print('[prvkey-check] failed!')

		# double-check pubkey
		if 1:
			#pubkey = str(bytes_to_hex(PublicKey.from_valid_secret(int_to_bytes_padded(prvkey)).format(flag_compress)))
			pubkey = getPubkey(prvkey,flag_compress)

			if pubkey != pubkey0:
				print('[pubkey#%s] %s' % (pow2bits,pubkey))
				print('[origin#%s] %s' % (pow2bits,pubkey0))
				print('[pubkey-check] failed!')

		# location in keyspace on the strip
		if 1:
			if flag_debug > -1:
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
		
		# finish stat
		printstr = '[i] %s j/s; %sj of %sj %.1f%%; DP T%s+W%s=%s+%s=%s; dp/kgr=%.1f' % (
			 prefSI(runjump/runtime)
			, runjump if runjump<10**3 else prefSI(runjump)
			, 2*Wsqrt if 2*Wsqrt < 10**3 else prefSI(2*Wsqrt)
			, (1.0*runjump/(2*Wsqrt))*100
			, HTmax,HWmax, lenT,lenW, lenT+lenW, 1.0*(lenT+lenW)/(HTmax+HWmax)
			)
		printstr += '' if (HTmax==1 and HWmax==1) else '; %s/%s repairs' % (nTrepair,nWrepair)
		printstr += '  '
		print(printstr)
		#print('[runtime]%s' % time_format(runtime))
		print('[runtime]%s' % time_format(runtime, (1,1,1,1,1,1,0,0)))


	print("[################################################]")

	#avgTime = (time.time()-starttime)/Ntimeit
	avgTime = 1.0 * sum(runtime for runtime in list_runtime) / len(list_runtime)
	avgJump = 1.0 * sum(runjump for runjump in list_runjump) / len(list_runjump)
	avgDPkg = 1.0 * sum(rundpkg for rundpkg in list_dpkgr) / len(list_dpkgr)
	#D = sum((xi - avgJump) ** 2 for xi in list_runjump)*1.0 / len(list_runjump)

	if Ntimeit > 1:
		##print('[(avg)jump] %.0f' % (avgJump) )
		#print('[(avg)jump] %s ' % (int(avgJump) if avgJump<10**3 else prefSI(avgJump)) )
		##print('[(avg)jum2] %.1f +/- %.1f' % (avgJump, (D/(len(list_runjump)-1))**0.5) )
		#print('[(avg)dpkg] %s ' % (int(avgDPkg) if avgDPkg<10**3 else prefSI(avgDPkg)) )
		#print('[(avg)time]%s' % time_format(avgTime, (1,1,1,1,1,1,1,0)) )

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
						,time_format( avgTime * xi , (1,1,1,1,1,1,1,0)) 
						,time_format( (avgTime * xi)/(avgJump/(2*Wsqrt)) , (1,1,1,1,1,1,1,0)) 
					) 
			)
		print("----------------------------------------------------------------------------------------------")
	else:
		pass

	print("[################################################]")
	print('[date] {}'.format(time.ctime()))
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
						,time_format( avgTime * xi , (1,1,1,1,1,1,1,0)) 
						,time_format( (avgTime * xi)/(avgJump/(2*Wsqrt)) , (1,1,1,1,1,1,1,0)) 
					) 
			)
		for i in xrange(pow2W,bitMax+1):
			xi = ((2**i)**0.5) / Wsqrt
			print('%s2^%03d |  %06s/ %06s |%032s /%032s |' % 
					(	'>' if i==pow2W else ' '
						,i
						,int(avgJump*xi) if int(avgJump*xi)<10**3 else prefSI(avgJump*xi)
						,int(2*(2**i)**0.5) if int(2*(2**i)**0.5)<10**3 else prefSI(2*(2**i)**0.5)
						,time_format( avgTime * xi , (1,1,1,1,1,1,1,0)) 
						,time_format( (avgTime * xi)/(avgJump/(2*Wsqrt)) , (1,1,1,1,1,1,1,0)) 
					) 
			)
		print("----------------------------------------------------------------------------------------------")

	print("[################################################]")
	#print('[date] {}'.format(time.ctime()))
	print('[exit] exit')
	#print('');raw_input('Press ENTER to continue...');print('')
	exit(0)
