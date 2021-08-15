module PositionTypes

__precompile__(true)

export DTE, IV, Options, Underlying, Position
export to_τ, get_ivol, set_DTE, set_ivol

mutable struct DTE{T<:Real}
    days::T
end

mutable struct IV{T<:Real}
    imp_vol::T
end

struct Options{T<:Real}
    isCall::Bool
    strike::T
    expiry::DTE
    ivol::IV
end

struct Underlying{T<:Real}
    price::T
end

struct Position{T<:Real, S}
    holding::Union{Options{T}, Underlying{T}}
    quantity::S
end

to_τ(x::Options) = x.expiry.days / 365.0 # We use ACT/365 here
get_ivol(x::Options) = x.ivol.imp_vol

function set_DTE(x::Options, new_days)
    @assert new_days != 0 "For Greek calculation DTE cannot be 0"
    x.expiry.days = new_days
end

function set_ivol(x::Options, new_ivol)
    @assert new_ivol != 0 "For Greek calculation IVol cannot be 0"
    x.ivol.imp_vol = new_ivol
end

end