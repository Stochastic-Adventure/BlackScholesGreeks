module GreekAnalyticalInv

__precompile__(true)

using Distributions

export blackScholesInvValue, blackScholesInvDelta, blackScholesInvGamma, blackScholesInvVega, blackScholesInvTheta
export blackScholesInvVanna, blackScholesInvVolga, blackScholesInvCharm
export blackScholesInvRho, blackScholesInvRho2

std_normal = Normal(0, 1)

## y = funding rate, corresponds to Espen's book, is q. 
## y = r - b, b = r - y

function blackScholesInvValue(isCall::Bool, spot::T, strike::T, τ::T, r::T, y::T, σ::T) where {T}
    sign = isCall ? 1 : -1

    d1 = log(spot / strike) / σ / sqrt(τ) + (r - y + 0.5 * σ) * sqrt(τ)
    d2 = d1 - σ * sqrt(τ)
    d3 = d2 - σ * sqrt(τ)

    sign * (exp(-r * τ) * cdf(std_normal, sign * d2) - exp((y - r + σ^2) * τ) * strike / spot * cdf(std_normal, sign * d3))
end

function blackScholesInvDelta(isCall::Bool, spot::T, strike::T, τ::T, r::T, y::T, σ::T) where {T}
    sign = isCall ? 1 : -1
    d1 = log(spot / strike) / σ / sqrt(τ) + (r - y + 0.5 * σ) * sqrt(τ)
    d2 = d1 - σ * sqrt(τ)
    d3 = d2 - σ * sqrt(τ)

    (exp(-r * τ) - exp((y - r) * (1 - σ) * τ)) / spot / σ / sqrt(τ) * pdf(std_normal, d2) + sign * exp((y - r + σ*σ) * τ) * strike / spot / spot * cdf(std_normal, sign * d3)
end

function blackScholesInvGamma(isCall::Bool, spot::T, strike::T, τ::T, r::T, y::T, σ::T) where {T}
    sign = isCall ? 1 : -1
    d1 = log(spot / strike) / σ / sqrt(τ) + (r - y + 0.5 * σ) * sqrt(τ)
    d2 = d1 - σ * sqrt(τ)
    d3 = d2 - σ * sqrt(τ)

    pdf_coeff1 = (1.5 + (r - y) / σ + log(spot / strike) / σ / σ / τ) * exp((y - r) * (1 - σ) * τ)
    pdf_coeff2 = (0.5 + (r - y) / σ + log(spot / strike) / σ / σ / τ) * exp(-r * τ)

    (pdf_coeff1 - pdf_coeff2) / spot / spot / σ / sqrt(τ) * pdf(std_normal, d2) - sign * exp((y - r + σ * σ) * τ) / (spot^3) * strike * 2 * cdf(std_normal, sign * d3)
end

function blackScholesInvVega(isCall::Bool, spot::T, strike::T, τ::T, r::T, y::T, σ::T) where {T}
    sign = isCall ? 1 : -1
    d1 = log(spot / strike) / σ / sqrt(τ) + (r - y + 0.5 * σ) * sqrt(τ)
    d2 = d1 - σ * sqrt(τ)
    d3 = d2 - σ * sqrt(τ)

    pdf_coeff1 = exp((y - r) * (1 - σ) * τ) * (1.5 * sqrt(τ) + log(spot / strike) / σ / σ / sqrt(τ))
    pdf_coeff2 = exp(-r * τ) * (0.5 * sqrt(τ) + log(spot / strike) / σ / σ / sqrt(τ))

    (pdf_coeff1 - pdf_coeff2) * pdf(std_normal, d2) - sign * 2 * exp((y - r + σ * σ) * τ) * strike / spot * σ * τ * cdf(std_normal, sign * d3)
end

function blackScholesInvTheta(isCall::Bool, spot::T, strike::T, τ::T, r::T, y::T, σ::T) where {T}
    sign = isCall ? 1 : -1
    d1 = log(spot / strike) / σ / sqrt(τ) + (r - y + 0.5 * σ) * sqrt(τ)
    d2 = d1 - σ * sqrt(τ)
    d3 = d2 - σ * sqrt(τ)

    pdf_coeff1 = exp((y - r) * (1 - σ) * τ) * ((y - r + 1.5 * σ) / 2 / sqrt(τ) + log(spot / strike) / 2 / σ / τ / sqrt(τ))
    pdf_coeff2 = exp(-r * τ) * ((y - r + 0.5 * σ) / 2 / sqrt(τ) + log(spot / strike) / 2 / σ / τ / sqrt(τ))

    sign * exp(-r * τ) * r * cdf(std_normal, sign * d2) + sign * exp((y - r + σ * σ) * τ) * strike / spot * (y - r + σ * σ) * cdf(std_normal, sign * d3) - (pdf_coeff1 - pdf_coeff2) * pdf(std_normal, d2)
end

function blackScholesInvVanna(isCall::Bool, spot::T, strike::T, τ::T, r::T, y::T, σ::T) where {T}
    sign = isCall ? 1 : -1
    d1 = log(spot / strike) / σ / sqrt(τ) + (r - y + 0.5 * σ) * sqrt(τ)
    d2 = d1 - σ * sqrt(τ)
    d3 = d2 - σ * sqrt(τ)

    pdf_coeff1 = exp((y - r) * (1 - σ) * τ) * (y - r) * sqrt(τ) / spot / σ
    pdf_coeff2 = (exp(-r * τ) -  exp((y - r) * (1 - σ) * τ)) / spot / σ / σ / sqrt(τ)
    pdf_coeff3 = exp((y - r) * (1 - σ) * τ) / spot * (1.5 * sqrt(τ) + log(spot / strike) / σ / σ / sqrt(τ))
    pdf_coeff4 = σ * pdf_coeff2 * (0.5 * sqrt(τ) + log(spot / strike) / σ / σ / sqrt(τ)) * ((r - y - 0.5 * σ) * sqrt(τ) + log(spot / strike) / σ / sqrt(τ))

    sign * 2 * strike * σ * τ / spot / spot * exp((y - r + σ * σ) * τ) * cdf(std_normal, sign * d3) + (pdf_coeff1 - pdf_coeff2 - pdf_coeff3 + pdf_coeff4) * pdf(std_normal, d2)
end

function blackScholesInvVolga(isCall::Bool, spot::T, strike::T, τ::T, r::T, y::T, σ::T) where {T}
    sign = isCall ? 1 : -1
    d1 = log(spot / strike) / σ / sqrt(τ) + (r - y + 0.5 * σ) * sqrt(τ)
    d2 = d1 - σ * sqrt(τ)
    d3 = d2 - σ * sqrt(τ)

    pdf_coeff1 = 2 * exp((y - r) * (1 - σ) * τ) * σ * τ * (1.5 * sqrt(τ) + log(spot / strike) / σ / σ / sqrt(τ))
    pdf_coeff2 = (0.5 * sqrt(τ) + log(spot / strike) / σ / σ / sqrt(τ)) * ((r - y - 0.5 * σ) * sqrt(τ) + log(spot / strike) / σ / sqrt(τ)) * (exp((y - r) * (1 - σ) * τ) * (1.5 * sqrt(τ) + log(spot / strike) / σ / σ / sqrt(τ)) - exp(-r * τ) * (0.5 * sqrt(τ) + log(spot / strike) / σ / σ / sqrt(τ)))
    pdf_coeff3 = 2 * log(spot / strike) / σ / σ / σ / sqrt(τ) * (exp(-r * τ) -  exp((y - r) * (1 - σ) * τ)) - exp((y - r) * (1 - σ) * τ) * (y - r) * τ * (1.5 * sqrt(τ) + log(spot / strike) / σ / σ / sqrt(τ))

    -2 * sign * exp((y - r + σ * σ) * τ) * strike * τ / spot * (1 + 2 * σ * σ * τ) * cdf(std_normal, sign * d3) + (pdf_coeff1 + pdf_coeff2 + pdf_coeff3) * pdf(std_normal, d2)
end

function blackScholesInvCharm(isCall::Bool, spot::T, strike::T, τ::T, r::T, y::T, σ::T) where {T}
    sign = isCall ? 1 : -1
    d1 = log(spot / strike) / σ / sqrt(τ) + (r - y + 0.5 * σ) * sqrt(τ)
    d2 = d1 - σ * sqrt(τ)
    d3 = d2 - σ * sqrt(τ)

    pdf_coeff1 = (exp(-r * τ) -  exp((y - r) * (1 - σ) * τ)) / 2 / spot / σ / τ / sqrt(τ)
    pdf_coeff2 = (r * exp(-r * τ) + (y - r) * (1 - σ) * exp((y - r) * (1 - σ) * τ)) / spot / σ / sqrt(τ)
    pdf_coeff3 = ((y - r + 1.5 * σ) / 2 / sqrt(τ) + log(spot / strike) / 2 / σ / τ / sqrt(τ)) * exp((y - r) * (1 - σ) * τ) / spot
    pdf_coeff4 = (exp(-r * τ) -  exp((y - r) * (1 - σ) * τ)) * ((y - r + 0.5 * σ) / 2 / sqrt(τ) + log(spot / strike) / 2 / σ / τ / sqrt(τ)) * ((y - r + 0.5 * σ) * sqrt(τ) - log(spot / strike) / σ / sqrt(τ)) / spot / σ / sqrt(τ)

    (pdf_coeff1 + pdf_coeff2 + pdf_coeff3 + pdf_coeff4) * pdf(std_normal, d2) - sign * exp((y - r + σ * σ) * τ) / spot / spot * strike * (y - r + σ * σ) * cdf(std_normal, sign * d3)
end

function blackScholesInvRho2(isCall::Bool, spot::T, strike::T, τ::T, r::T, y::T, σ::T) where {T}
    sign = isCall ? 1 : -1
    d1 = log(spot / strike) / σ / sqrt(τ) + (r - y + 0.5 * σ) * sqrt(τ)
    d2 = d1 - σ * sqrt(τ)
    d3 = d2 - σ * sqrt(τ)

    (exp((y - r) * (1 - σ) * τ) - exp(-r * τ)) * sqrt(τ) * pdf(std_normal, d2) - sign * exp((y - r + σ * σ) * τ) * strike * τ / spot * cdf(std_normal, sign * d3)
end

function blackScholesInvRho(isCall::Bool, spot::T, strike::T, τ::T, r::T, y::T, σ::T) where {T}
    sign = isCall ? 1 : -1
    d1 = log(spot / strike) / σ / sqrt(τ) + (r - y + 0.5 * σ) * sqrt(τ)
    d2 = d1 - σ * sqrt(τ)
    d3 = d2 - σ * sqrt(τ)

    sign * strike / spot * τ * exp((y - r + σ * σ) * τ) * cdf(std_normal, sign * d3) - sign * exp(-r * τ) * τ * cdf(std_normal, sign * d2) + (exp(-r * τ) - exp((y - r) * (1 - σ) * τ)) * sqrt(τ) * pdf(std_normal, d2)
end

end