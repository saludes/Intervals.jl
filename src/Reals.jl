#module Reals

using RadixIntervals, Intervals, RadixStreams

immutable Real
	r::Function
end

typealias Q Rational{BigInt}


function Real(a::Approx)
	function approx(e)
		upto!(a, ifloor(log(a, e)))
		convert(Interval{Q}, a)
	end # approx
	Real(approx)
end

function Real(r::RadixNumber)
	function approx(e)
		a = Approx(r)
		upto!(a, ifloor(log(a, e)))
		convert(Interval{Q}, a)
	end # approx
	Real(approx)
end

function Real{T<:Number}(q::Rational{T})
	function rational(_)
		Interval(big(q))
	end # rational
	Real(rational)
end

function check(r::Real)
	function ck()
		e = 1.0e0
		i0 = nothing
		while true
			i = r.r(e)
			produce((e, diam(i) <= 2e))
			e /= 2.0
		end
	end
	return Task(ck)
end



# Test
using Sqrt

n2 = RadixNumber(2, DEC)
b = sqrt(Approx(n2))
s2 = Real(n2)

#end # module