module RadixIntervals

using Intervals, RadixStreams


import Base.convert

function convert{T<:Number}(::Type{Interval{T}}, a::Approx)
	q = convert(T, a)
	eps = 1//big(a.radix)^(-exponent(a))
	q >= 0 ? Interval(q, q + eps) : Interval(q - eps, q)
end


end # module

