module alphaPredictableComponents

using LinearAlgebra, Printf, Statistics
using DataFrames, CoupledFields, FFTW, StatsBase

export modelf, prca, pca, CVfn
export Rsq1, Rsq2, Rsq3
export Rsq, omega, iota
export LCCA
export CoDa

include("CoDa.jl")


"""
    parzen(tau::OrdinalRange{Int64}, M::Int64)

Parzen filter
"""
function parzen(tau::OrdinalRange{Int64}, M::Int64)
    # lags =-50:50; f=  parzen(lags, 50)
    # p=plot(x=lags, y=f, Geom.line); draw(PNG(6.6inch, 3.3inch), p)

    τ = abs.(tau)
    i = τ .>= M/2   
    z =   (i).*(2(1 .- τ./M).^3) +  (.!i).*(1 .- 6(τ./M).^2 + 6(τ./M).^3)        
    z *= 1.0/sum(z)
    return z
end

"""
    modelf(par::Vector, X::Matrix, modelfn::Function)

Find the components of ``X`` using `modelfn` and parameters `par`.
par[1] = ``α``, par[2] = ``τ``
"""
function modelf(par::Vector{Float64}, X::Matrix{Float64}, modelfn::Function)
    Z0 = CoDa.alphaTransform(X, par[1])
    Z = Z0 .- mean(Z0,dims=1)
    lagmax = floor(Int, size(Z,1)/8)
    return modelfn(par, Z, lagmax)
end

"""
    prca(par::Vector, X::Matrix, lagmax::Int)

Predictable Components Analysis.
"""
function prca(par::Vector{Float64}, X::Matrix{Float64}, lagmax::Int)
    local n::Int64
    local p::Int64
    local Cxx_fac_inv::Matrix{Float64}
    
    n,p = size(X)
    Cxx_fac_inv = cov(X)^(-0.5)
    Zx = X * Cxx_fac_inv
#    lagmax = floor(Int, n/4)
    filt = parzen(0:lagmax, lagmax) 
    mc = StatsBase.crosscov(Zx, Zx, 0:lagmax, demean=false)
    M = zeros(Float64, p, p)
    for k in 1:lagmax
        m1 = mc[1+k,:,:]
        M += filt[1+k] * m1 * m1'    
    end
    
    d,U = eigen(Symmetric(M))
    W = Cxx_fac_inv * reverse(U,dims=2)
    Rx = X * W
    A1 = zeros(Float64, p, 4)
    for j in 1:4
        A1[:,j] = LCCA(Rx[:,j:j], Zx, Int(par[2]))
    end
    A = Cxx_fac_inv * A1
    Ry = X * A
    return ModelObj(W, Rx, A, Ry, reverse(d), par, "αPrC")
end

"""
    pca(par::Vector, X::Matrix, lagmax::Int)

Principal Components Analysis
"""
function pca(par::Vector{Float64}, X::Matrix{Float64}, lagmax::Int)
    
    Csc = diagm(0=>vec(1 ./ std(X, dims=1)))
    Z = X * Csc 
    F, d, V = svd(Z)
#    W = Csc * V
    Rx = F./std(F,dims=1)
    W = (X'X) \ (X'Rx[:,1:4])  
    Ry = X*W
    return ModelObj(W, Rx, W, Ry, d.^2, par, "αPC")
end

"""
    CVfn(parm::Matrix, X::Matrix, modelfn::Function, cvmetric::Function)

Cross-validation
"""
function CVfn(parm::Matrix{Float64}, X::Matrix{Float64}, modelfn::Function, cvmetric::Function) 
    
    # CV metric parameters
    # par[2] = τ lags
    
    # Model parameters
    # par[1] = alpha

    nrg = size(parm, 1)
    n, p = size(X)
     
    # Setup
    ta = 1:floor(Int, n/2)
    tb = setdiff(1:n,ta)

    lagmax = floor(Int, n/8)
    
    CVmetric = zeros(Float64, nrg)
    for i in 1:nrg
        percent = round(100*i/nrg)
        @printf " %03d%%" percent

        par = parm[i,:]
        τ = Int(par[2])
        Z = CoDa.alphaTransform(X, par[1])
        zm1 = mean(Z[ta,:], dims=1); # zm2 = mean(Z[tb,:],1)
    
        Z1train = Z[ta,:].-zm1 
        Z1test = Z[tb,:].-zm1 

        model1 = modelfn(par, Z1train, lagmax)
#        modelfn==pca && println("αPCA")
#        modelfn==pca && (model1 = Afn(model1, Z1train, τ))
        CVmetric[i] = cvmetric(model1, Z1test, τ)
        println([parm[i:i,:] CVmetric[i] ])
    end
    parm = [parm CVmetric]
    imax = argmax(parm[:,end])
    return parm[imax,:]

end


"""
    Rsq1(model::ModelObj, Xtest::Matrix, τ::Int)

Predictive R²
"""
function Rsq1(model::ModelObj, Xtest::Matrix{Float64}, τ::Int)
    Rcv = Xtest * model.W[:,1]
    Tcv = Xtest * model.A[:,1]
    return cor(Rcv[1:(end-τ),1], Tcv[(1+τ):end,1])^2
end

"""
    Rsq2(model::ModelObj, Xtest::Matrix, τ::Int)

Explanation of variance R²
"""
function Rsq2(model::ModelObj, Xtest::Matrix{Float64}, τ::Int)
    Rcv = Xtest * model.W[:,1]
    numer = (Xtest'Rcv) / (Rcv'Rcv) * (Rcv'Xtest) 
    denom = Xtest'Xtest
    return tr(numer)/tr(denom)
end    

"""
    Rsq3(model::ModelObj, Xtest::Matrix, τ::Int)


"""
function Rsq3(model::ModelObj, Xtest::Matrix{Float64}, τ::Int)
    Rcv = Xtest * model.W[:,1]
    numer  = sum(abs2, Rcv / (Rcv'Rcv) * (Rcv'Xtest))
    denom = sum(abs2, Xtest)
    return numer/denom
end    


function Afn(model::ModelObj, X::Matrix{Float64}, τ::Int)
    
    p = size(X,2)
    Cxx_fac_inv = cov(X)^(-0.5)
    Zx = X * Cxx_fac_inv
    Rx = X * model.W
    A1 = zeros(Float64, p, 4)
    for j in 1:4
        A1[:,j] = LCCA(Rx[:,j:j], Zx, τ)
    end
    A = Cxx_fac_inv * A1
    Ry = X * A
    return ModelObj(model.W, Rx, A, Ry, model.evals, model.pars, "αPCA")
end


"""
    LCCA(X::Matrix, Y::Matrix, τ::Int)

Linear Canonical Correlation Analysis between ``X`` and ``Y`` at lead time ``τ``
"""
function LCCA(X::T, Y::T, τ::Int) where T<:Matrix{Float64}
    nr = size(X,1)-τ
    Xa = X[1:nr,:]
    Ya = Y[(1+τ):end,:]
    M = cor(Xa, Ya)
    U, d, V = svd(M)
#    sv = cor(Xa*U[:,1], Ya*V[:,1])
    return V[:,1]
end


"""
    biplot(model::ModelObj, jv::Vector{Int}, Y::Matrix, label::Vector{Symbol})

Biplot 
"""
function biplot(model::ModelObj, jv::Vector{Int}, Y::Matrix{Float64}, label::Vector{Symbol})

    Z = CoDa.alpha(Y, model.pars[1])
    R = model.R[:,jv]    
    B = cov(R, zscore(Z,1))
    
    j = string(model.method,jv[2])
    return DataFrame(label=string.(label), B1=B[1,:], B2=B[2,:], x0=0.0, y0=0.0, a=string("α=",model.pars[1]), j=j)
end    


"""
    Rsq(R::Matrix, X0::Matrix)

Explanation of variance R²
"""
function Rsq(R::Matrix{Float64}, X0::Matrix{Float64})
    numer = mapslices(Rj->tr((X0'Rj)/(Rj'Rj)*(Rj'X0)), R, dims=1)
    denom = tr(X0'X0)
    return vec(numer)./denom
end    
    

"""
    omega(X::Matrix)

Spectral entropy
"""
function omega(X::Matrix{Float64})
    Xfft = fft(X,1)[2:end,:]
    sdf = real(Xfft .* conj(Xfft))
    n = size(Xfft,1)
    sdf[sdf.<0.0] .= 0.0
    sdf ./= sum(sdf,dims=1)
    o =  sdf .* log.(n, sdf)
    o[isinf.(o)] .= 0.0
    ιfca = 1.0 .+ sum(o,dims=1)
    return vec(ιfca)
end

omega(X::Matrix{Float64}, X0::Matrix{Float64}) = omega(X)


"""
    iota(model::ModelObj, iotafn::Function, jv::Vector, X0::Matrix, pos::Vector)

Iota
"""
function iota(model::ModelObj, iotafn::Function, jv::Vector{Int}, X0::Matrix{Float64}, pos::Vector)
    R = model.R[:,jv]
    ι = round.(Int, 100 .*iotafn(R, X0))
    istr = iotafn==omega ? string.("Ω = ",ι, "%") : string.(ι, "%")
    j = model.method .*string.(jv)  
    return DataFrame(j=j, label=istr, x=pos[1], y=pos[2])
end    



end # module
