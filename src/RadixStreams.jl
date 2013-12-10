module RadixStreams


export
	RadixNumber

export
	Packet,
	Approx, upto!, exponent, toInt,
	explode, implode, split, pop!, hmerge

export
	BIN, DEC, OCT, HEX

typealias Radix Uint8
typealias Exponent Int32
typealias Packet (Exponent, BigInt)

maxexp = typemax(Exponent)

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


# Extract digits upto! an exponent

function pop!(a::Approx)
	p = consume(a.stream)
	push!(a.packets, p)
	return p
end

function upto!(a::Approx, n::Integer)
	if n >= exponent(a)
		return # Do nothing
	end
	for (e,m) = a.stream
		push!(a.packets, @pkt e m)
			if e <= n
				break
			end
	end
	if istaskdone(a.stream) && n < exponent(a) # exact: add zeros
		push!(a.packets, @pkt n 0)
	end
end




# Truncation

function trunc(n::RadixNumber, e::Exponent)
 	a = Approx(n)
 	upto!(a, n)
	convert(RadixNumber, a)
end

import Base.convert


function convert(::Type{RadixNumber}, a::Approx)
	function truncated()
		for p=a.packets
			produce(p)
		end
	end # truncated
	RadixNumber(a.radix, truncated)
end


# Printing

int2char(r::Uint8) = char(r > 9 ? int('a') + r -10 : int('0') + r)
char2int(c::Char) = c >= 'a' ? int(c) - int('a') + 10 : int(c) - int('0') 

function fromInt(n::Integer,b::Radix)
	digs = Radix[]
	while n > 0
		n,r = divmod(n, b)
		#s = "$(int2char(r))$s"
		push!(digs, r)
	end
	return digs
end



function toInt(s::Array{Integer,1},b::Radix)
	v = 0
	for d = s
		v = v*b + d # char2int(d)
	end
	return v
end

toInt(s::String, b::Radix) = toInt(map(char2int, convert(Array{Char,1},s)))


function explode(a::Approx; forever=false)
	function exploded()
		exp = maxexp
		for (e,m) = a.stream
			digs = Packet[]
			for p = zip(e+[0:1000], fromInt(m ,a.radix))
				push!(digs, p)
			end
			while maxexp > exp > length(digs) + e
				produce((exp, 0))
				exp -= 1
			end
			for p = reverse(digs)
				produce(p)
				exp = p[1]
			end
		end
		while forever # Continue giving 0's
			exp -= 1
			produce(@pkt exp 0)
		end
	end # exploded
	Approx(a.radix, Task(exploded), Packet[])
end

function implode(l::Array{Packet,1}, b::Radix)
	if isempty(l)
		return (maxexp,0)
	end
	value = big(0)
	e1 = minimum(map(p -> p[1], l))
	for (e,d) = l
		value += d * b^(e - e1)
	end
	return (e1,value)
end

function hmerge(a::Approx, n)
	pks = splice!(a.packets, 1:n)
	unshift!(a.packets, implode(pks, a.radix))
	return a.packets[1]
end

function split(a::Approx, n::Integer)
	function splat()
		fp = Packet[]
		for p = explode(a, forever=true).stream
			push!(fp, p)
			if p[1] % n == 0
				produce(implode(fp, a.radix))
				empty!(fp)
			end
		end
	end
	Approx(a.radix, Task(splat), Packet[])
end



import Base.print

function print(io::IO, x::RadixNumber)
	a = Approx(x)
	b = a.radix
	exponent = maxexp
	
	function format(digs,e1)
		e2 = e1 + length(digs)
		if e1 < 0 <= e2
			"$(digs[1:e2]).$(digs[e2+1:end])"
		else
			digs
		end
	end

	pr(s) = print(io, s)
	function lpadto(s, n)
		if n <= length(s)
			return s
		end
		zeros = "0"^(n-length(s))
		return "$zeros$s"
	end

	pr("$b) ")
	for (e1,n) = a.stream
		ds = fromInt(n, b)
		dgs = isempty(ds) ? "0" : convert(String, map(int2char, ds))
		if exponent < maxexp
			dgs = lpadto(dgs, exponent - e1) # add significant zeros
		end
		pr(format(dgs, e1))
		exponent = e1
	end
end

function print(io::IO, a::Approx)
	function pr(s)
			prev ? print(io, " + ") : nothing
			prev = true
			print(io, s)
	end
	prev = false
	for (e,m) = a.packets
		if e == 0
			pr(m == 0 ? "" : "$m")
		else
			pr("$m*$(a.radix)^$e")
		end
	end
	istaskdone(a.stream) ? println(io, "") : pr("...\n")
end


function convert{T<:Number}(::Type{Rational{T}}, a::Approx)
	num = zero(T)
	b = convert(T, a.radix)
	e0 = exponent(a)
	if e0 == maxexp # Not started
		return nothing
	end
	for (e,m) = a.packets
		num += m * b^(e-e0)
	end
	e0 >=0 ? Rational(num * b^e0) : Rational(num, b^(-e0))
end


end # module