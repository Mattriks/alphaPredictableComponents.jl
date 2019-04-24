
module CoDa

export alpha, alphaTransform, Helmert
export ilr, clr 

"""
    Helmert(n)

Calculate ``n×n`` Helmert matrix
"""
function Helmert(n)
    H = zeros(n,n)
  for i = 1:n, j = 1:n
      if i==1
        H[i,j] = 1.0/sqrt(n)
      elseif j < i
        H[i,j] = 1.0/sqrt(i*(i-1))
      elseif i==j
        H[i,j] = (1.0-i)/sqrt(i*(i-1))
      end
  end
    return -H
end

"""
    alphaTransform(X::Matrix, α::Float64)

Apply alpha-transform with Helmert product to ``X``
"""
function alphaTransform(X::Matrix{Float64}, α::Float64)
    D = size(X,2)
    H = Helmert(D)[2:end,:]
    Xᵃ = X.^α
    Xs = sum(Xᵃ, dims=2)    
    U = (D*Xᵃ./Xs .- 1)./α
    return U*transpose(H)    
end


"""
    alpha(X::Matrix, α::Float64)

Apply alpha-transform without Helmert product
"""
function alpha(X::Matrix{Float64}, α::Float64)
    D = size(X,2)
#    H = Helmert(D)[2:end,:]
    Xᵃ = X.^α
    Xs = sum(Xᵃ,dims=2)    
    U = (D*Xᵃ./Xs .- 1)./α
    return U    
end


"""
    ilr(X::Matrix)

Isometric logratio transform
"""
function ilr(X::Matrix{Float64})
    D = size(X,2)  
    H = Helmert(D)[2:end,:]
    return log.(X) * transpose(H)  
end    


"""
    clr(X::Matrix)

Centered logratio transform
"""
function clr(X::Matrix{Float64})
    D = size(X,2)
    ratio = X ./ (prod(X,2).^(1/D))
    return log.(ratio)
end


end