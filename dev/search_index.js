var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#alphaPredictableComponents.jl-1",
    "page": "Home",
    "title": "alphaPredictableComponents.jl",
    "category": "section",
    "text": ""
},

{
    "location": "#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": ""
},

{
    "location": "#alpha-Transform-1",
    "page": "Home",
    "title": "alpha Transform",
    "category": "section",
    "text": ""
},

{
    "location": "#Predictability-1",
    "page": "Home",
    "title": "Predictability",
    "category": "section",
    "text": "iota_Pr=fracWleft(sum_tau=1^infty Sigma_YXSigma_X^-1Sigma_XY right) WW Sigma_XW..."
},

{
    "location": "#Guide-1",
    "page": "Home",
    "title": "Guide",
    "category": "section",
    "text": "Pages = joinpath.(\"man\", [\"Example1.md\"])\nDepth = 1"
},

{
    "location": "#Library-1",
    "page": "Home",
    "title": "Library",
    "category": "section",
    "text": "Pages = [joinpath(\"lib\",\"library.md\")]"
},

{
    "location": "#Documentation-1",
    "page": "Home",
    "title": "Documentation",
    "category": "section",
    "text": "The documentation was built using Documenter.jl.println(\"Documentation built with Julia $(VERSION).\") # hide"
},

{
    "location": "man/Example1/#",
    "page": "Guide",
    "title": "Guide",
    "category": "page",
    "text": ""
},

{
    "location": "man/Example1/#Example-1-1",
    "page": "Guide",
    "title": "Example 1",
    "category": "section",
    "text": ""
},

{
    "location": "lib/library/#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "lib/library/#Library-1",
    "page": "Library",
    "title": "Library",
    "category": "section",
    "text": "Pages = [\"library.md\"]"
},

{
    "location": "lib/library/#Types-1",
    "page": "Library",
    "title": "Types",
    "category": "section",
    "text": "ModelObj"
},

{
    "location": "lib/library/#alphaPredictableComponents.CVfn-Tuple{Array{Float64,2},Array{Float64,2},Function,Function}",
    "page": "Library",
    "title": "alphaPredictableComponents.CVfn",
    "category": "method",
    "text": "CVfn(parm::Matrix, X::Matrix, modelfn::Function, cvmetric::Function)\n\nCross-validation\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.LCCA-Union{Tuple{T}, Tuple{T,T,Int64}} where T<:Array{Float64,2}",
    "page": "Library",
    "title": "alphaPredictableComponents.LCCA",
    "category": "method",
    "text": "LCCA(X::Matrix, Y::Matrix, τ::Int)\n\nLinear Canonical Correlation Analysis between X and Y at lead time τ\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.Rsq-Tuple{Array{Float64,2},Array{Float64,2}}",
    "page": "Library",
    "title": "alphaPredictableComponents.Rsq",
    "category": "method",
    "text": "Rsq(R::Matrix, X0::Matrix)\n\nExplanation of variance R²\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.Rsq1-Tuple{CoupledFields.ModelObj,Array{Float64,2},Int64}",
    "page": "Library",
    "title": "alphaPredictableComponents.Rsq1",
    "category": "method",
    "text": "Rsq1(model::ModelObj, Xtest::Matrix, τ::Int)\n\nPredictive R²\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.Rsq2-Tuple{CoupledFields.ModelObj,Array{Float64,2},Int64}",
    "page": "Library",
    "title": "alphaPredictableComponents.Rsq2",
    "category": "method",
    "text": "Rsq2(model::ModelObj, Xtest::Matrix, τ::Int)\n\nExplanation of variance R²\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.Rsq3-Tuple{CoupledFields.ModelObj,Array{Float64,2},Int64}",
    "page": "Library",
    "title": "alphaPredictableComponents.Rsq3",
    "category": "method",
    "text": "Rsq3(model::ModelObj, Xtest::Matrix, τ::Int)\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.iota-Tuple{CoupledFields.ModelObj,Function,Array{Int64,1},Array{Float64,2},Array{T,1} where T}",
    "page": "Library",
    "title": "alphaPredictableComponents.iota",
    "category": "method",
    "text": "iota(model::ModelObj, iotafn::Function, jv::Vector, X0::Matrix, pos::Vector)\n\nIota\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.modelf-Tuple{Array{Float64,1},Array{Float64,2},Function}",
    "page": "Library",
    "title": "alphaPredictableComponents.modelf",
    "category": "method",
    "text": "modelf(par::Vector, X::Matrix, modelfn::Function)\n\nFind the components of X using modelfn and parameters par. par[1] = α, par[2] = τ\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.omega-Tuple{Array{Float64,2}}",
    "page": "Library",
    "title": "alphaPredictableComponents.omega",
    "category": "method",
    "text": "omega(X::Matrix)\n\nSpectral entropy\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.pca-Tuple{Array{Float64,1},Array{Float64,2},Int64}",
    "page": "Library",
    "title": "alphaPredictableComponents.pca",
    "category": "method",
    "text": "pca(par::Vector, X::Matrix, lagmax::Int)\n\nPrincipal Components Analysis\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.prca-Tuple{Array{Float64,1},Array{Float64,2},Int64}",
    "page": "Library",
    "title": "alphaPredictableComponents.prca",
    "category": "method",
    "text": "prca(par::Vector, X::Matrix, lagmax::Int)\n\nPredictable Components Analysis.\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.biplot-Tuple{CoupledFields.ModelObj,Array{Int64,1},Array{Float64,2},Array{Symbol,1}}",
    "page": "Library",
    "title": "alphaPredictableComponents.biplot",
    "category": "method",
    "text": "biplot(model::ModelObj, jv::Vector{Int}, Y::Matrix, label::Vector{Symbol})\n\nBiplot \n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.parzen-Tuple{OrdinalRange{Int64,S} where S,Int64}",
    "page": "Library",
    "title": "alphaPredictableComponents.parzen",
    "category": "method",
    "text": "parzen(tau::OrdinalRange{Int64}, M::Int64)\n\nParzen filter\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.CoDa.Helmert-Tuple{Any}",
    "page": "Library",
    "title": "alphaPredictableComponents.CoDa.Helmert",
    "category": "method",
    "text": "Helmert(n)\n\nCalculate nn Helmert matrix\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.CoDa.alpha-Tuple{Array{Float64,2},Float64}",
    "page": "Library",
    "title": "alphaPredictableComponents.CoDa.alpha",
    "category": "method",
    "text": "alpha(X::Matrix, α::Float64)\n\nApply alpha-transform without Helmert product\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.CoDa.alphaTransform-Tuple{Array{Float64,2},Float64}",
    "page": "Library",
    "title": "alphaPredictableComponents.CoDa.alphaTransform",
    "category": "method",
    "text": "alphaTransform(X::Matrix, α::Float64)\n\nApply alpha-transform with Helmert product to X\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.CoDa.clr-Tuple{Array{Float64,2}}",
    "page": "Library",
    "title": "alphaPredictableComponents.CoDa.clr",
    "category": "method",
    "text": "clr(X::Matrix)\n\nCentered logratio transform\n\n\n\n\n\n"
},

{
    "location": "lib/library/#alphaPredictableComponents.CoDa.ilr-Tuple{Array{Float64,2}}",
    "page": "Library",
    "title": "alphaPredictableComponents.CoDa.ilr",
    "category": "method",
    "text": "ilr(X::Matrix)\n\nIsometric logratio transform\n\n\n\n\n\n"
},

{
    "location": "lib/library/#Functions-1",
    "page": "Library",
    "title": "Functions",
    "category": "section",
    "text": "Modules = [alphaPredictableComponents, alphaPredictableComponents.CoDa]\nOrder   = [:function]"
},

]}
