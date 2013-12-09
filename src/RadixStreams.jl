module RadixStreams


export
	RadixNumber,
	Approx, upto, exponent, 
	BIN, DEC, OCT, HEX

typealias Radix Uint8
typealias Exponent Int32
typealias Packet (Exponent, BigInt)

maxexp= typemax(Exponent)

# Usual radices
for (i,n) in [(:BIN,2), (:DEC, 10), (:HEX, 16), (:OCT, 8)]
	@eval ($i = uint8($n))
end




# Like Python's divmod
divmod(n, m) = (div(n,m), rem(n,m))


immutable RadixNumber
	radix::Radix
	digits::Function
end

type Approx
	radix::Radix
	stream::Task
	packets::Array{Packet,1}
end

Approx(n::RadixNumber) = Approx(n.radix, Task(n.digits), Packet[])

exponent(a::Approx) = isempty(a.packets) ? maxexp : a.packets[end][1]

macro pkt(e,n)
	:( (int32($e), big($n)) )
end



getpkt(a::Approx) = begin
	(e,n) = consume(a.stream)
	p = @pkt e n
	append!(a.packets, p)
	return p
end


macro ckradix(r)
	:(1 < $r < typemax(Radix) ? uint8($r) : error("Invalid radix: ", $(string(r))))
end

function RadixNumber(x::Integer, r::Integer)
	b = @ckradix r
	function intdigits()
		produce(@pkt 0 x)
		return
	end
	RadixNumber(b, intdigits)
end

function RadixNumber(x::Rational, r::Integer)
	b = @ckradix r
	m = den(x)

	function rdivide(e,n)
		q,r = divmod(n, m)
		@assert q < b
		produce(@pkt e q)
		if r == 0
			return
		else
			rdivide(e-1, b*r)
		end
	end

	function rational()
		ip,rest = divmod(num(x), m)
		produce((int32(0),big(ip)))
		if rest == 0
			return
		end
		rdivide(-1, b*rest)
	end
	RadixNumber(b, rational)
end


# Extract digits upto an exponent
function upto(a::Approx, n::Integer)
	if n < exponent(a)
		for (e,m) = a.stream
			append!(a.packets, [@pkt e m])
			if e <= n
				break
			end
		end
	end
end


# Truncation

function trunc(n::RadixNumber, e::Exponent)
 	a= Approx(n)
 	upto(a, n)
	function truncated()
		for p=a.packets
			produce(p)
		end
 	end
 	RadixNumber(n.radix, truncated)
end


# Printing

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

function print(io::IO, x::RadixNumber)
	a = Approx(x)
	b = a.radix
	exponent = maxexp
	print(io, "$b) ")
	for (e1,n) in a.stream
		dgs = fromInt(n, b)
		e2 = e1 + length(dgs)
		if exponent < maxexp
			print(io, convert(String, repeat(['0'], inner=[exponent-e2])))
		end
		if e1 < 0 <= e2
			print(io, "$(dgs[1:e2]).$(dgs[e2+1:end])")
		else
			print(io, dgs)
		end
		exponent = e1
	end
end

function print(io::IO, a::Approx)
	for (e,m) = a.packets
		if e == 0
			print(io, m == 0 ? "" : "$m + ")
		else
			print(io, "$m*$(a.radix)^$e  + ")
		end
	end
	println("...")
end

import Base.convert

function convert(Rational, a::Approx)
	num = zero(BigInt)
	b = big(a.radix)
	e0 = exponent(a)
	for (e,m) = a.packets
		num += m * b^(e-e0)
	end
	e0 >=0 ? Rational(num * b^e0) : Rational(num, b^(-e0))
end

end # module