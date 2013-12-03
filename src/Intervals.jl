# This file is part of Intervals.jl.
#
# Intervals.jl is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# Intervals.jl is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Intervals.jl; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# Also add information on how to contact you by electronic and paper mail.
#
# Copyright (C) 2013 Alessandro Andrioni
# Copyright (C) 2013 Jordi Saludes


module Intervals

export
    Interval,
    QQ,
    # bisect,
    # blow,
    # diam,
    # diam_abs,
    # diam_rel,
    # mag,
    mid
    # mig,
    # isbounded

import
    Base: precision, string, print, show, showcompact, promote_rule,
        promote, convert, +, *, -, /, exp, isinf, isnan, nan, inf, sqrt,
        square, exp, exp2, expm1, cosh, sinh, tanh, sech, csch, coth, inv,
        sqrt, cbrt, abs, log, log2, log10, log1p, sin, cos, tan, sec,
        csc, acos, asin, atan, acosh, asinh, atanh, isempty, union,
        intersect, in, cmp, inv

typealias QQ Rational{BigInt}
typealias IntervalTypes Union(QQ)
typealias SmallFloat Union(Float32, Float64)

immutable Interval{T<:IntervalTypes}
    left::T
    right::T

    function Interval(left::T, right::T)
        if left > right
            error("Invalid: $left > $right")
        end
        new(left, right)
    end
end

Interval{T<:IntervalTypes}(x::T) = Interval{T}(x, x)
Interval{T<:IntervalTypes}(left::T, right::T) = Interval{T}(left, right)
Interval(l::Rational, r::Rational) = Interval(big(l), big(r))

# Conversion and promotion related functions
# Conversions to Interval
convert(::Type{Interval{QQ}}, x::Number) = Interval(convert(QQ, x))

# Conversions from Interval

#convert(::Type{QQ}, i::Interval) = convert(QQ, mid(i))
#for to in (Int8, Int16, Int32, Int64, Uint8, Uint16, Uint32, Uint64, BigInt, Float32)
#    @eval begin
#        function convert(::Type{$to}, x::Interval)
#            convert($to, convert(QQ, x))
#        end
#    end
#end

promote_rule{T<:Number,S<:Number}(::Type{T}, ::Type{Interval{S}}) = Interval{promote_type(T,S)}


# Basic operations
function +(x::Interval, y::Interval)
    l = x.left + y.left
    r = x.right + y.right
    Interval(l, r)
end
function -(x::Interval, y::Interval)
    l = x.left - y.right
    r = x.right - y.left
    Interval(l, r)
end
function *(x::Interval, y::Interval)
    if x.left >= 0
        if y.left >= 0
            l = *(x.left, y.left)
            r = *(x.right, y.right)
        elseif y.right <= 0
            l = *(x.right, y.left)
            r = *(x.left, y.right)
        else
            l = *(x.right, y.left)
            r = *(x.right, y.right)
        end
    elseif x.right <= 0
        if y.left >= 0
            l = *(x.left, y.right)
            r = *(x.right, y.left)
        elseif y.right <= 0
            l = *(x.right, y.right)
            r = *(x.left, y.left)
        else
            l = *(x.left, y.right)
            r = *(x.left, y.left)
        end
    else
        if y.left >= 0
            l = *(x.left, y.right)
            r = *(x.right, y.right)
        elseif y.right <= 0
            l = *(x.right, y.left)
            r = *(x.left, y.left)
        else
            l = min(*(x.left, y.right), *(x.right, y.left))
            r = max(*(x.left, y.left), *(x.right, y.right))
        end
    end
    Interval(l, r)
end

*(x::Number, i::Interval) = (@show promote_type(typeof(x),typeof(i)); *(promote(x,i)...))
# TODO: Reimplement directly without using inv
function /(x::Interval, y::Interval)
    i = inv(y)
    x * i
end

# General functions
function mid(x::Interval)
    m = x.left + x.right
    if isinf(m)
        return (x.left / 2) + (x.right / 2) # ??
    else
        return m / 2
    end
end

function inv{T<:IntervalTypes}(x::Interval{T})
    l = one(T)/x.right
    r = one(T)/x.left
    Interval(l, r)
end

# Printing-related functions
string(x::Interval) = "[$(string(x.left)), $(string(x.right))]"
print(io::IO, x::Interval) = print(io, string(x))
show(io::IO, x::Interval) = print(io, string(x))
showcompact(io::IO, x::Interval) = print(io, string(x))


end # module
