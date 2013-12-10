module Sqrt

	using RadixStreams, Intervals, RadixIntervals

import Base.sqrt

function sqrt(a::Approx)
	B = a.radix
	a2 = split(a, 2)
	sq = zero(BigInt)
	function squared()
		(e,m) = pop!(a2)
		while m > 0
			res = map(r -> (r, m - r * (2 * a.radix * sq + r)), [0:B])
			(d,m) = filter(x -> x[2]>=0, res)[end]
			@assert d < B
			a2.packets[end] = (e,m)
			np = (div(e,2), d)
			sq = sq * B + d
			produce(np)
			pop!(a2)
			(e,m) = hmerge(a2, 2)
		end
	end # squared
	Approx(a.radix, Task(squared), Packet[])
end


# Test sqrt(2)
typealias Q Rational{BigInt}
exp = 50

n2 = RadixNumber(2, HEX)
b = sqrt(Approx(n2))
upto!(b, -exp)
i = convert(Interval{Q}, b)
@assert crossing(x -> x^2 - 2, i)
@assert diam(i) < 2epsilon(b)

end # module