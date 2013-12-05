#module RadixRing

#export
#	Radix, divide,
#	BIN, DEC, OCT, HEX

typealias Radix Uint8

for (i,n) in [(:BIN,2), (:DEC, 10), (:HEX, 16), (:OCT, 8)]
	@eval ($i = uint8($n))
end

immutable RadixStream
	radix::Radix
	digits::Function
end

macro ckradix(r)
	:(1 < $r < 26 ? uint8($r) : error("Invalid radix: ", $(string(r))))
end

function RadixStream(x::Integer, r::Integer)
	b = @ckradix r
	function f()
		produce(x)
		return
	end
	RadixStream(b, f)
end

function RadixStream(x::Rational, r::Integer)
	b = @ckradix r
	function f()
		divide(num(x), den(x), b)
	end
	RadixStream(b, f)
end


divmod(n, m) = (div(n,m), rem(n,m))


function divide(n::Integer, m::Integer, b::Radix)
	function divi(n)
		if n == 0
			return
		end
		q,r = divmod(n, m)
		produce(q)
		divi(b*r)
	end
	divi(n)
end


shift(x::Radix,n::Integer=1...) = Radix(x.radix, x.mantisa*x.radix^n, uint16(x.exponent+n))


int2char(r::Uint8) = char(r > 9 ? int('a') + r -10 : int('0') + r)
char2int(c::Char) = c >= 'a' ? int(c) - int('a') + 10 : int(c) - int('0') 

function fromInt(n::Integer,b::Uint8)
	if n==0
		return "0"
	end
	s = ""
	while n>0
		n,r = div(n,b), uint8(rem(n,b))
		s = "$(int2char(r))$s"
	end
	return s
end

function toInt(s::String,b::Uint8)
	v = 0
	for d=s
		v = v*b + char2int(d)
	end
	return v
end

import Base.print

function print(io::IO, x::RadixStream)
	ds = Task(x.digits)
	fi(n) = fromInt(n, x.radix)
	print(io, fi(consume(ds)))
	d = consume(ds)
	if d == nothing # Integer
		return
	end
	print(io,".")
	print(io, fi(d))
	for d=ds
		print(io, fi(d))
	end
end

#end # module