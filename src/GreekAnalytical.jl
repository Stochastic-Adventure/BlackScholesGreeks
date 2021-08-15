module GreekAnalytical

__precompile__(true)

using Distributions

export blackScholesValue, blackScholesDelta, blackScholesGamma, blackScholesVega, blackScholesTheta
export blackScholesVanna, blackScholesVolga, blackScholesSpeed, blackScholesCharm, blackScholesColor, blackScholesZomma
export blackScholesRho, blackScholesRhoFutures, blackScholesTouchProb, blackScholesDriftlessTheta

std_normal = Normal(0, 1)

function blackScholesValue(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    sign = isCall ? 1 : -1

    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    d2 = d1 - σ * sqrt(τ)

    sign * (spot * exp((b-r) * τ) * cdf(std_normal, sign * d1) - strike * exp(-r * τ) * cdf(std_normal, sign * d2))
end

function blackScholesDelta(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    if isCall
        exp((b-r) * τ) * cdf(std_normal, d1)
    else
        exp((b-r) * τ) * (cdf(std_normal, d1) - 1.0)
    end
end

function blackScholesGamma(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    pdf(std_normal, d1) * exp((b-r) * τ) / spot / σ / sqrt(τ)
end

function blackScholesVega(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    spot * exp((b-r) * τ) * pdf(std_normal, d1) * sqrt(τ)
end

function blackScholesTheta(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    sign = isCall ? 1 : -1

    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    d2 = d1 - σ * sqrt(τ)

    term1 = -spot * exp((b-r) * τ) * pdf(std_normal, d1) * σ / 2 / sqrt(τ)
    term2 = -sign * (b-r) * spot * exp((b-r) * τ) * cdf(std_normal, sign * d1)
    term3 = -sign * r * strike * exp(-r * τ) * cdf(std_normal, sign * d2)

    term1 + term2 + term3
end

function blackScholesVanna(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    d2 = d1 - σ * sqrt(τ)

    -d2 * pdf(std_normal, d1) / σ * exp((b-r) * τ)
end

function blackScholesCharm(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    sign = isCall ? 1 : -1

    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    d2 = d1 - σ * sqrt(τ)

    term1 = pdf(std_normal, d1) * (b / σ / sqrt(τ) - d2 / 2 / τ)
    term2 = sign * (b - r) * cdf(std_normal, sign * d1)

    -exp((b-r) * τ) * (term1 + term2)
end

function blackScholesVolga(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    d2 = d1 - σ * sqrt(τ)

    blackScholesVega(isCall, spot, strike, τ, r, b, σ) * d1 * d2 / σ
end

function blackScholesSpeed(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))

    -1.0 * blackScholesGamma(isCall, spot, strike, τ, r, b, σ) / spot * (1 + d1 / σ / sqrt(τ))
end

function blackScholesZomma(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    d2 = d1 - σ * sqrt(τ)

    blackScholesGamma(isCall, spot, strike, τ, r, b, σ) * (d1 * d2 - 1.0) / σ
end

function blackScholesColor(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    d2 = d1 - σ * sqrt(τ)

    blackScholesGamma(isCall, spot, strike, τ, r, b, σ) * (r - b + b * d1 / σ / sqrt(τ) + (1 - d1 * d2) / 2 / τ)
end

function blackScholesDriftlessTheta(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    println(pdf(std_normal, d1))
    -spot * pdf(std_normal, d1) * σ / 2 / sqrt(τ)
end

function blackScholesRho(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    sign = isCall ? 1 : -1

    d1 = (log(spot / strike) + (b + 0.5 * σ^2) * τ) / (σ * sqrt(τ))
    d2 = d1 - σ * sqrt(τ)

    sign * τ * strike * exp(-r * τ) * cdf(std_normal, sign * d2)
end

function blackScholesRhoFutures(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where{T}
    -1.0 * τ * blackScholesValue(isCall, spot, strike, τ, r, b, σ)
end

function blackScholesTouchProb(isCall::Bool, spot::T, strike::T, τ::T, r::T, b::T, σ::T) where {T}
    if isCall
        if strike <= spot
            return 1.0
        end
    else
        if strike >= spot
            return 1.0
        end
    end

    sign = isCall ? -1 : 1

    μ = (b - σ^2 / 2) / (σ * σ)
    λ = sqrt(μ^2 + 2 * r / (σ * σ))
    z = log(strike / spot) / σ / sqrt(τ) + λ * σ * sqrt(τ)

    (strike / spot) ^ (μ + λ) * cdf(std_normal, sign * z) + (strike / spot) ^ (μ - λ) * cdf(std_normal, sign * (z - 2 * λ * σ * sqrt(τ)))
end

end