module SpreadBuilder

__precompile__(true)

using ..GreekAnalytical
using ..PositionTypes: Position
export portfolio_greeks, portfolio_PnL

function portfolio_greeks(portf::AbstractVector{Position}, underlying_price::T, r::T, b::T) where {T<:Real}
    ref_maturity = minimum(x->x.DTE, portf) / 365.0

    portf_delta = 0.0
    portf_gamma = 0.0
    portf_vega = 0.0
    portf_theta = 0.0
    portf_vanna = 0.0
    portf_volga = 0.0
    portf_charm = 0.0
    portf_color = 0.0

    for pos in portf
        if pos.holding isa Options
            option = pos.holding
            quantity = pos.quantity
            portf_delta += blackScholesDelta(option.isCall, underlying_price, option.strike, to_τ(option), r, b, get_ivol(option)) * quantity
            portf_gamma += blackScholesGamma(option.isCall, underlying_price, option.strike, to_τ(option), r, b, get_ivol(option)) * quantity
            portf_theta += blackScholesTheta(option.isCall, underlying_price, option.strike, to_τ(option), r, b, get_ivol(option)) * quantity
            portf_vega += sqrt(to_τ(option) / ref_maturity) * blackScholesVega(option.isCall, underlying_price, option.strike, to_τ(option), r, b, get_ivol(option)) * quantity
            portf_vanna += sqrt(to_τ(option) / ref_maturity) * blackScholesVanna(option.isCall, underlying_price, option.strike, to_τ(option), r, b, get_ivol(option)) * quantity
            portf_volga += sqrt(to_τ(option) / ref_maturity) * blackScholesVolga(option.isCall, underlying_price, option.strike, to_τ(option), r, b, get_ivol(option)) * quantity
            portf_charm += blackScholesCharm(option.isCall, underlying_price, option.strike, to_τ(option), r, b, get_ivol(option)) * quantity
            portf_color += blackScholesColor(option.isCall, underlying_price, option.strike, to_τ(option), r, b, get_ivol(option)) * quantity
        else
            quantity = pos.quantity
            portf_delta += pos.quantity
        end
    end
end

function portfolio_PnL(portf::AbstractVector{Position}, initial_price::T, underlying_price::T, r::T, b::T) where {T<:Real}
    pnl = 0.0

    for pos in portf
        if pos.holding isa Options
            option = pos.holding
            quantity = pos.quantity
            pnl += (blackScholesValue(option.isCall, underlying_price, option.strike, to_τ(option), r, b, get_ivol(option)) - blackScholesValue(option.isCall, initial_price, option.strike, to_τ(option), r, b, get_ivol(option))) * quantity
        else
            pnl += underlying_price - initial_price
        end
    end
end

end