
import itertools
import math as m 

class TimeStampedVector():
	
	#Format is timestamp first, then scalar as [[time,scalar],[time,scalar],[time,scalar]]
	
	id_iter = itertools.count()
	
	def get_self_name(self):
		return self.Name
	
	def get_identity_from_index(self, index):
		return self.Data[index]
	
	def get_timestamp_from_index(self, index):
		return self.Data[index][0]
	
	def get_scalar_from_index(self, index):
		return self.Data[index][1]
	
	def set_self_name(self,string_name):
		self.id = next(self.id_iter)
		self.Name = "%s #%s" % (str(string_name), self.id)
	
	def set_identity_from_index(self, index, identity=[0,0]):
		self.Data[index] = identity
	
	def set_timestamp_from_index(self, index, value=0):
		self.Data[index][0] = value
		
	def set_scalar_from_index(self, index, value=0):
		self.Data[index][1] = value
	
	def len(self):
		return len(self.Data)
	
	#Identity is a short name for the list with the time and scalar, as in [time,scalar] to append to the Data
	def append_identity(self, identity):
		self.Data.append(identity)
	
	def create_sin_wave(self, start_time=0, end_time=20, resolution=10, freq=3, magnitude=1):
		q = TimeStampedVector(None,"create_sin_wave")
		for i in range(m.floor(end_time-start_time)*resolution):
			time = start_time + i/(end_time*resolution)
			q.append_identity([time,m.sin(time)])
		return q
	
	def create_square_wave(self, start_time=0, end_time=20, resolution=10, freq=3, magnitude=1):
		q = TimeStampedVector(None,"create_square_wave_q")
		for i in range(m.floor((end_time-start_time)*resolution)):
			time = start_time + i/(end_time-start_time*resolution)
			scalar = -magnitude
			if time % (1/freq) < (1/(2*freq)):
				scalar = -scalar
			q.append_identity([time,scalar])
		return q
	
	def add_waves(self, TimeStampedVectorA, TimeStampedVectorB):
		q = TimeStampedVector(None,"add_waves")
		for i, e in enumerate(TimeStampedVectorA.Data):
			q.append_identity(	[TimeStampedVectorA.get_timestamp_from_index(i),
								TimeStampedVectorA.get_scalar_from_index(i)
								+TimeStampedVectorB.get_scalar_from_index(i)]
								)
		return q
	
	#Print Data visually to the terminal as formatted "%s: %s|" % (time, number of spaces according to scalar) if padding is True
	def print_to_terminal(self, scalar=1, offset=0, padding=True):
		for i in self.Data:
			if padding:
				t = format(m.floor(i[0]*200)/200, "#.3g") + ": "
			else:
				t = str(i[0]) + ": "
			for j in range(m.floor(i[1]*scalar + offset)):
				t += " "
			t += "|"
			print(t)
	
	def __init__(self, predefined_two_d_vector=None, name=None):
		
		if isinstance(predefined_two_d_vector,type(None)):
			self.Data = [[0,0],[0,0],[0,0]]
		else:
			self.Data = predefined_two_d_vector
		
		self.id = next(self.id_iter)
		
		if isinstance(name,type(None)):
			self.Name = "Undefined #%s" % (self.id)
		else:
			self.Name = "%s #%s" % (name,self.id)


class Complex():
	
	def add(self, ComplexA, ComplexB):
		
		q = Complex(ComplexA.Real + ComplexB.Real, ComplexA.Imaginary + ComplexB.Imaginary)
		
		return q

	def sub(self, ComplexA, ComplexB):
		
		q = Complex(ComplexA.Real - ComplexB.Real, ComplexA.Imaginary - ComplexB.Imaginary)
		
		return q
	
	def mul(self, ComplexA, ComplexB):
		
		q = Complex(ComplexA.Real * ComplexB.Real - ComplexA.Imaginary * ComplexB.Imaginary, ComplexA.Real * ComplexB.Imaginary + ComplexA.Imaginary * ComplexB.Real)
		
		return q
		
	def div(self, ComplexA, ComplexB):
		
		q = Complex((ComplexA.Real * ComplexB.Real + ComplexA.Imaginary * ComplexB.Imaginary)/(ComplexB.Real * ComplexB.Real + ComplexB.Imaginary * ComplexB.Imaginary), 
					(ComplexA.Imaginary * ComplexB.Real - ComplexA.Real * ComplexB.Imaginary)/(ComplexB.Real * ComplexB.Real + ComplexB.Imaginary * ComplexB.Imaginary))
		
		return q
	
	def distance_to_origin(self, ComplexN):
		
		return m.sqrt(ComplexN.Real**2+ComplexN.Imaginary**2)
	
	def __init__(self, real=0, imaginary=0):
		
		self.Real = real
		self.Imaginary = imaginary


def dft(TimeStampedVectorN,start_freq=0,end_freq=20,resolution=5):
	#Discrete Fourier Tranformation implemented by Kim C. Chemnitz
	#resolution refers to how many frequencies/iterations to compute per unit of time
	
	complex = Complex()
	
	q = TimeStampedVector(None,"dft_q")
	
	freq = start_freq
	
	origin_time = TimeStampedVectorN.get_timestamp_from_index(0)
	
	length = TimeStampedVectorN.len()
	
	while freq < end_freq:
		
		t = Complex()
		
		for i, e in enumerate(TimeStampedVectorN.Data):
			
			time = TimeStampedVectorN.get_timestamp_from_index(i) - origin_time
			
			angle = 2*m.pi*time*freq
			scalar = TimeStampedVectorN.get_scalar_from_index(i)
			
			temp = Complex(scalar*m.sin(angle), scalar*m.cos(angle))
			
			t = complex.add(t, temp)
		
		q.append_identity([freq,complex.distance_to_origin(t)])
		
		freq += 1/resolution
	
	return q
