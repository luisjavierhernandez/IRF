module Irf

greet() = print("Hello World!")

import Colors
#import PyPlot
import Polynomials

export rgb2arrayBasinsUnivariateRationalFunctions, polynomialToNewtonRationalFunction, fixedIndeterminatePointsofaPair


function sphereBijection_pairofcomplexes(
    twotuple::Tuple{Complex{Float64},Complex{Float64}})
    z=twotuple[1]
    t=twotuple[2]
    den=real(conj(t)*t + conj(z)*z)
        if den==0.0
            point=[0.0,0.0,0.0]
        else
            point=[real((conj(z)*t + conj(t)*z)/den),
            real((1im*(conj(z)*t - conj(t)*z))/den),
            real((-conj(t)*t + conj(z)*z)/den)]
        end
    return point
end

#Example:
sphereBijection_pairofcomplexes((1.0+0.0*im, 0.0+1.0*im))

function chordalMetric(twotuple::Tuple{Complex{Float64},Complex{Float64}},
    twotuple1::Tuple{Complex{Float64},Complex{Float64}})
    norma=vector->sqrt(vector[1]^2+vector[2]^2+vector[3]^2)
    return norma((sphereBijection_pairofcomplexes(twotuple)-
    sphereBijection_pairofcomplexes(twotuple1)))
end

#Example
chordalMetric((1.0+0.0*im, 0.0+1.0*im),(0.0+0.0*im, 1.0+0.0*im))

#Using two parametrization we can cover all the projective plane



function rectangle(xinterval::Tuple{Float64,Float64}=(-1.5,1.5),
                    yinterval::Tuple{Float64,Float64}=(-1.5,1.5),
                    expprecision=10)
    tol=1.0/(2^expprecision)
    a=xinterval[1]
    b=xinterval[2]
    c=yinterval[1]
    d=yinterval[2]
    red=[(complex(r,i),complex(1.0)) for i=d:-tol:c, r=a:tol:b]
    return red
end

rectanglecompila=rectangle((-2.56,2.56),(-2.,2.),3)


function invertedrectangle(xinterval::Tuple{Float64,Float64}=(-1.5,1.5),
                    yinterval::Tuple{Float64,Float64}=(-1.5,1.5),
                    expprecision=10)
    tol=1.0/(2^expprecision)
    a=xinterval[1]
    b=xinterval[2]
    c=yinterval[1]
    d=yinterval[2]
    red=[(complex(1.0), complex(r,i)) for i=d:-tol:c, r=a:tol:b]
    return red
end

invertedrectanglecompila=invertedrectangle((-2.56,2.56),(-2.56,2.56),3)



function homogeneousNormalization_pairofcomplexexNullity(twotuple::Tuple{T,S},
    nullitytolerance::Int64=15) where {T <: Number, S <: Number}
    tt1=complex(twotuple[1])
    tt2=complex(twotuple[2])
    den=abs(tt1)+abs(tt2)
    if den<1.0/10^nullitytolerance
        hpoint=(complex(0.0),complex(0.0))
    else
    hpoint=(tt1/den, tt2/den)
    return hpoint
    end
end

#Examples
homogeneousNormalization_pairofcomplexexNullity((1.0+0.0*im, 0.0+1.0*im),9)

homogeneousNormalization_pairofcomplexexNullity((0.000000000000000000001+0.0*im,
        0.0+0.0*im),9)


function samenumberofcoeficcients(list1,list2)
          le1=length(list1)
          le2=length(list2)
          if le1==le2
            pairoflists=[complex(list1),complex(list2)]
          end
          if le1<le2
            len1=le1
            lista1=list1
            while len1<le2
              lista1= push!(lista1,complex(0.0))
              len1=len1+1
            end
              pairoflists=[lista1,list2]
          end
          if le1>le2
            len2=le2
            lista2=list2
            while len2<le1
              lista1= push!(lista2,complex(0.0))
              len2=len2+1
            end
              pairoflists=[complex(list1),complex(lista2)]
          end
          le=max(le1,le2)
          firstlist=pairoflists[1]
          secondlist=pairoflists[2]
          while firstlist[le]==complex(0.0) && secondlist[le]==complex(0.0) && le>1
              firstlist=deleteat!(firstlist,le)
              secondlist=deleteat!(secondlist,le)
              le=le-1
            end
            return [firstlist, secondlist]
end

function normalize_pairhpolynomials(numerator,
                denominator)
                norma=sum(abs,numerator)+sum(abs,denominator)
                newnum=numerator*(1.0/norma)
                newden=denominator*(1.0/norma)
            return [newnum,newden]
end

function polynomialToNewtonRationalFunction(coefficienlist)
            pol=Polynomials.Poly(coefficienlist)
            polx=Polynomials.Poly([0.0,1.0])
            numerator=Polynomials.polyder(pol)*polx-pol
            denominator=Polynomials.polyder(pol)
            coeffnumerator=complex([numerator[i] for i in 0:1:length(numerator)])
            coeffdenominator=complex([denominator[i] for i in 0:1:length(denominator)])
            nundem=samenumberofcoeficcients(coeffnumerator, coeffdenominator)
            nundemnorm=normalize_pairhpolynomials(nundem[1],nundem[2])
            return  nundemnorm
end

rfcompila=polynomialToNewtonRationalFunction([-1.0,0.0,0.0,1.0])
coefficientlistnumcompila=rfcompila[1]
coefficientlistdencompila=rfcompila[2]

function bivariatepolyfunction(coefficientlist,
     u, t, d::Int)
    coefficientlistcomplex=complex(coefficientlist)
    ff=coefficientlistcomplex[1]*t^d
    for i in 2:length(coefficientlist)
        ff=ff+coefficientlist[i]*u^(i-1)*t^(d-i+1)
    +i
    end
    return ff
end

function bivariatepolyfunction(coefficientlist::Array{T,1},
      u::S, t::R, d::Int) where {T<:Number,S<:Number,R<:Number}
     coefficientlistcomplex=complex(coefficientlist)
     ff=coefficientlistcomplex[1]*t^d
     for i in 2:length(coefficientlist)
         ff=ff+coefficientlist[i]*u^(i-1)*t^(d-i+1)
     +i
     end
     return ff
end




function pairoffunctions(coefficientlistnum,
        coefficientlistden)
        ln=length(coefficientlistnum)-1
        ld=length(coefficientlistden)-1
        d=max(ln,ld)
        fff(u, t)=
        bivariatepolyfunction(coefficientlistnum,u,t,d)
        ggg(u, t)=
        bivariatepolyfunction(coefficientlistden,u,t,d)
        return fff, ggg
end



function pairoffunctions(coefficientlistnum::Array{A,1},
         coefficientlistden::Array{B,1}) where {A<:Number, B<:Number}
         ln=length(coefficientlistnum)-1
         ld=length(coefficientlistden)-1
         d=max(ln,ld)
         function fff(u::S, t::R) where {S<:Number,R<:Number}
         bivariatepolyfunction(coefficientlistnum,u,t,d)
            end
         function ggg(u::S, t::R) where {S<:Number,R<:Number}
         bivariatepolyfunction(coefficientlistden,u,t,d)
            end
        return fff, ggg
end



function rationalFunction(hpair::Tuple{Function,Function},
        twotuple::Tuple{Complex{Float64},Complex{Float64}},
        nullitytolerance::Int64=15)
        c=twotuple[1]
        d=twotuple[2]
        F=hpair[1]
        G=hpair[2]
        cnew=F(c,d)
        dnew=G(c,d)
        hresult=homogeneousNormalization_pairofcomplexexNullity((cnew,dnew),
        nullitytolerance)
    return  hresult
end




hpaircompila=pairoffunctions(coefficientlistnumcompila,
    coefficientlistdencompila)

rationalFunction(hpaircompila,(1.0+0.0*im, 0.0+1.0*im),9)





## The thoery of the map theta and fixed and indeterminate points.

function fixedIndeterminatePointsofaPair(coefficientlistnum,
        coefficientlistden)
        numden=samenumberofcoeficcients(coefficientlistnum,coefficientlistden)
        coefficientlistnumcomplex=numden[1]
        coefficientlistdencomplex=numden[2]
        ffgg=pairoffunctions(coefficientlistnumcomplex, coefficientlistdencomplex)
        norma=sum(abs,coefficientlistnumcomplex)+sum(abs,coefficientlistdencomplex)
        if norma==0.0
            newfixed=[(complex(0.0),complex(0.0))] ; morenewfixed=[];fixindpoints=["All the points are indetermination point points"]
            return newfixed, morenewfixed,  fixindpoints
        elseif length(coefficientlistnumcomplex)==1 && length(coefficientlistdencomplex)==1
            newfixed=[];morenewfixed=[];
                fixindpoints=[(coefficientlistnumcomplex[1],coefficientlistdencomplex[1])]
            return newfixed, morenewfixed, fixindpoints
        elseif abs(coefficientlistnumcomplex[1])==0.0 &&
                coefficientlistdencomplex[length(coefficientlistdencomplex)]==0.0 &&
                sum(i->abs(coefficientlistnumcomplex[i]-coefficientlistdencomplex[i-1]),
                2:length(coefficientlistdencomplex))==0.0
                newfixed=[(complex(0.0),complex(0.0))]
                morenewfixed=[(complex(1.0),complex(0.0))]
                fixindfixed=["All the points are fixed points"]
            return newfixed, morenewfixed, fixindfixed
        else
            if abs(ffgg[2](complex(1.0),complex(0.0)))==0.0 && abs(ffgg[1](complex(1.0),complex(0.0)))!=0.0
                newfixedpre=[(complex(0.0),complex(0.0))] ; morenewfixedpre=[(complex(1.0),complex(0.0))]; fixindfixedpre=[]
            elseif abs(ffgg[2](complex(1.0),complex(0.0)))==0.0 && abs(ffgg[1](complex(1.0),complex(0.0)))==0.0
                    newfixedpre=[(complex(0.0),complex(0.0))] ; morenewfixedpre=[]; fixindfixedpre=[(complex(1.0),complex(0.0))]
            end
            newfixed=newfixedpre; morenewfixed=morenewfixedpre; fixedindpoints=fixindfixedpre
            Polynumerator=Polynomials.Poly(coefficientlistnumcomplex)
            Polydenominator=Polynomials.Poly(coefficientlistdencomplex)
            Polyx=Polynomials.Poly([0.0+ 0.0* im,1.0+0.0*im])
            lookingzeros= Polynumerator- Polydenominator * Polyx
            fix=Polynomials.roots(lookingzeros)
            le=length(fix)
            mas=[(complex(fix[i]),1.0+0.0*im) for i in 1:le]
                if length(fixedindpoints)==0
                fixindfixed=mas
                else
                fixindfixed=prepend(fixedindpoints[1],mas)
                end
            return newfixedpre, morenewfixed, fixindfixed
        end
end

fixedIndeterminatePointsofaPairex=fixedIndeterminatePointsofaPair(coefficientlistnumcompila,coefficientlistdencompila)








#We can stop the process taking into account the chordal metric

function newstep(hpair::Tuple{Function,Function},
                    iter::Int64,
                    nullitytolerance::Int64,
                    cauchytolerance::Int64,
                    hpoint::Tuple{Complex{Float64},Complex{Float64}})
    point = hpoint
    tol=1.0/10^(cauchytolerance)
    number = 0
    imagepoint = rationalFunction(hpair,point,nullitytolerance)
    while (chordalMetric(point, imagepoint) > tol) && (number < iter)
    point = imagepoint
    imagepoint = rationalFunction(hpair,point,nullitytolerance)
    number=number+1
    end
    return [imagepoint,number]
end

newstep(hpaircompila, 25, 15, 3, (1.0+0.0*im, 0.0+1.0*im))

function newstep(hpair::Tuple{Function,Function},
        iter::Int64,
        nullitytolerance::Int64,
        cauchytolerance::Int64,
    hrectangle::Array{Tuple{Complex{Float64},Complex{Float64}},2})
        size1=size(hrectangle)
        result=[newstep(hpair, iter,nullitytolerance, cauchytolerance, hrectangle[i,j])
        for i=1:size1[1], j=1:size1[2]]
    return result
end

newstep(hpaircompila,25,15, 3, rectanglecompila)

function positionuptotolerance(fixedPointList::Array{Tuple{Complex{Float64},
                            Complex{Float64}},1},
                            nullitytolerance::Int64,
                            approachtolerance::Int64,
                            twotuple::Tuple{Complex{Float64},Complex{Float64}})
    pos=0
    it=1
    le = length(fixedPointList)
    tol=1/10^(approachtolerance)
    while (it < le+1)
        if (chordalMetric(twotuple, fixedPointList[it]) < tol)
            pos = it
        end
    it=it+1
    end
    return convert(Int64,pos)
end

function positionuptotolerance(fixedPointList::Array{Tuple{Complex{Float64},
                            Complex{Float64}},1},
                            nullitytolerance::Int64,
                            approachtolerance::Int64,
                            hrectangle::Array{Tuple{Complex{Float64},Complex{Float64}},2})
                            size1=size(hrectangle)
                            result=[positionuptotolerance(fixedPointList,  nullitytolerance, approachtolerance,hrectangle[i,j])
                            for i=1:1size1[1], j=1:1:size1[2]]
                        return result
end




positionuptotolerance(fixedIndeterminatePointsofaPairex[3],15,3,rectanglecompila)

function iteration_upto_tolerances(hpair::Tuple{Function,Function},
        iter::Int64,
        nullitytolerance::Int64,
        cauchytolerance::Int64,
        twotuple::Tuple{Complex{Float64},Complex{Float64}})
        endpoint_iterations=newstep(hpair, iter, nullitytolerance, cauchytolerance, twotuple)
        iterations=convert(Int64,endpoint_iterations[2])
    return iterations
end

iteration_upto_tolerances(hpaircompila,11, 3, 2, (1.0+0.0*im, 0.0+1.0*im))

function iteration_upto_tolerances(hpair::Tuple{Function,Function},
        iter::Int64,
        nullitytolerance::Int64,
        cauchytolerance::Int64,
        hrectangle::Array{Tuple{Complex{Float64},Complex{Float64}},2})
        size1=size(hrectangle)
        result=[iteration_upto_tolerances(hpair, iter, nullitytolerance, cauchytolerance,hrectangle[i,j])
        for i=1:1size1[1], j=1:1:size1[2]]
    return result
end

iteration_upto_tolerances(hpaircompila,11,15,3, rectanglecompila)



function position_iteration_upto_tolerances(hpair::Tuple{Function,Function},
        fixedPointList::Array{Tuple{Complex{Float64},Complex{Float64}},1},
        iter::Int64,
        tolerances::Tuple{Int64,Int64,Int64},
        twotuple::Tuple{Complex{Float64},Complex{Float64}})
    nullitytolerance=tolerances[1]
    cauchytolerance=tolerances[2]
    approachtolerance=tolerances[3]
    endpoint_iterations=newstep(hpair, iter, nullitytolerance, cauchytolerance, twotuple)
    endpoint=endpoint_iterations[1]
    iterations=convert(Int64,endpoint_iterations[2])
    pos=positionuptotolerance(fixedPointList, nullitytolerance, approachtolerance, endpoint)
    return pos, iterations
    endpoint
end

position_iteration_upto_tolerances(hpaircompila,fixedIndeterminatePointsofaPairex[3], 11, (15,3,2), (1.0+0.0*im, 0.0+1.0*im))

function position_iteration_upto_tolerances(hpair::Tuple{Function,Function},
        fixedPointList::Array{Tuple{Complex{Float64},Complex{Float64}},1},
        iter::Int64,
        tolerances::Tuple{Int64,Int64,Int64},
        hrectangle::Array{Tuple{Complex{Float64},Complex{Float64}},2})
    nullitytolerance=tolerances[1]
    cauchytolerance=tolerances[2]
    approachtolerance=tolerances[3]
    size1=size(hrectangle)
    result=[position_iteration_upto_tolerances(hpair,fixedPointList, iter,
        tolerances,
        hrectangle[i,j]) for i=1:size1[1], j=1:size1[2]]
    return result
end

position_iteration_upto_tolerances(hpaircompila,fixedIndeterminatePointsofaPairex[3],25, (15,3,2), rectanglecompila)



function fixedPointListex_matrixpositioninterations_RationalFunction(
        coefficientlistnum,
        coefficientlistden,
        expresolution::Int=8,
        iteration_max::Int=25,
        tolerances::Tuple{Int64,Int64,Int64}=(15,3,3),
        modeldomain::AbstractString="localrectangle",
        rectanglesidesdomain::Tuple{Float64,Float64,Float64,Float64}=(-1.5,1.5,-1.5,1.5))
    nullitytolerance=tolerances[1]
    cauchytolerance=tolerances[2]
    approachtolerance=tolerances[3]
    xinterv=(rectanglesidesdomain[1],rectanglesidesdomain[2])
    yinterv=(rectanglesidesdomain[3],rectanglesidesdomain[4])
    if modeldomain=="localrectangle"
        rect=rectangle(xinterv,  yinterv, expresolution)
    elseif modeldomain=="invertedlocalrectangle"
        rect=invertedrectangle(xinterv,  yinterv, expresolution)
    end
        hpair=pairoffunctions(coefficientlistnum, coefficientlistden)
        coefficientlistnumcomplex=complex(coefficientlistnum)
        coefficientlistdencomplex=complex(coefficientlistden)
    fixedPointListexcomplete=fixedIndeterminatePointsofaPair(coefficientlistnumcomplex,
        coefficientlistdencomplex)
    fixedPointListexcompletea=fixedPointListexcomplete[1]
    if length(fixedPointListexcomplete[2])==0
    fixedPointListexcompleteb=fixedPointListexcompletea
    else
    fixedPointListexcompleteb=union(fixedPointListexcompletea,
        [(complex(1.0),complex(0.0))])
    end
    fixedPointListexample=union(fixedPointListexcompleteb,fixedPointListexcomplete[3])
    positer=position_iteration_upto_tolerances(hpair,fixedPointListexample,
    iteration_max, tolerances, rect)
    return   fixedPointListexcomplete, positer
end

fixedPointListex_matrixpositioninterations_RationalFunction(coefficientlistnumcompila, coefficientlistdencompila, 3)


function auxfixedPointListex_matrixpositioninterations_RationalFunction(
        coefficientlistnum,
        coefficientlistden,
        expresolution::Int=8,
        iteration_max::Int=25,
        tolerances::Tuple{Int64,Int64,Int64}=(15,3,3),
        modeldomain::AbstractString="localrectangle",
        rectanglesidesdomain::Tuple{Float64,Float64,Float64,Float64}=(-1.5,1.5,-1.5,1.5))
    nullitytolerance=tolerances[1]
    cauchytolerance=tolerances[2]
    approachtolerance=tolerances[3]
    xinterv=(rectanglesidesdomain[1],rectanglesidesdomain[2])
    yinterv=(rectanglesidesdomain[3],rectanglesidesdomain[4])
    if modeldomain=="localrectangle"
        rect=rectangle(xinterv,  yinterv, expresolution)
    elseif modeldomain=="invertedlocalrectangle"
        rect=invertedrectangle(xinterv,  yinterv, expresolution)
    end
        hpair=pairoffunctions(coefficientlistnum, coefficientlistden)
        coefficientlistnumcomplex=complex(coefficientlistnum)
        coefficientlistdencomplex=complex(coefficientlistden)
    fixedPointListexcomplete=fixedIndeterminatePointsofaPair(coefficientlistnumcomplex,
        coefficientlistdencomplex)
    fixedPointListexcompletea=fixedPointListexcomplete[1]
    if length(fixedPointListexcomplete[2])==0
    fixedPointListexcompleteb=fixedPointListexcompletea
    else
    fixedPointListexcompleteb=union(fixedPointListexcompletea,
        [(complex(1.0),complex(0.0))])
    end
    fixedPointListexample=union(fixedPointListexcompleteb,fixedPointListexcomplete[3])
    positer=position_iteration_upto_tolerances(hpair,fixedPointListexample,
    iteration_max, tolerances, rect)
    aux=[length(fixedPointListexcomplete[1]),length(fixedPointListexcomplete[2]),length(fixedPointListexcomplete[3]),convert(Int64,iteration_max)]
    return   aux, positer
end



auxfixedPointListex_matrixpositioninterations_RationalFunction(rfcompila[1],rfcompila[2],8)

function auxfixedPointListex__RationalFunction(
        coefficientlistnum,
        coefficientlistden,
        iteration_max::Int=25)
        coefficientlistnumcomplex=complex(coefficientlistnum)
        coefficientlistdencomplex=complex(coefficientlistden)
    fixedPointListexcomplete=fixedIndeterminatePointsofaPair(coefficientlistnumcomplex,
        coefficientlistdencomplex)
    aux=[length(fixedPointListexcomplete[1]),length(fixedPointListexcomplete[2]),length(fixedPointListexcomplete[3]),convert(Int64,iteration_max)]
    return   aux
end

auxfixedPointListex__RationalFunction(rfcompila[1],rfcompila[2],25)

#the process can also be stoped using the function absolute value of numerator and denomonator

function image__rationalFunction(hpair,
        twotuple::Tuple{Complex{Float64},Complex{Float64}},
        nullitytolerance::Int64=15)
        c=twotuple[1]
        d=twotuple[2]
        F=hpair[1]
        G=hpair[2]
        cnew=F(c,d)
        dnew=G(c,d)
        hresult=(cnew,dnew)
        return   hresult
end



function abs_image__rationalFunction(hpair,
        twotuple::Tuple{Complex{Float64},Complex{Float64}},
        nullitytolerance::Int64=15)
        c=twotuple[1]
        d=twotuple[2]
        F=hpair[1]
        G=hpair[2]
        cnew=F(c,d)
        dnew=G(c,d)
        hresult=(cnew,dnew)
        ab=(abs(cnew)+abs(dnew))
        return   ab
end




function iterationFinalPoint__hpoint(hpair,
                    nullitytolerance::Int64,
                    firtsiternumber::Int64,
                    hpoint::Tuple{Complex{Float64},Complex{Float64}})
    tolnull=1.0/10^(nullitytolerance)
    zerozero=(complex(0.0),complex(0.0))
    number = 0
    norpoint=homogeneousNormalization_pairofcomplexexNullity(hpoint,
        nullitytolerance)
    imagepoint=image__rationalFunction(hpair,norpoint,nullitytolerance)
    norimagepoint=homogeneousNormalization_pairofcomplexexNullity(imagepoint,
        nullitytolerance)
    while number < firtsiternumber
    norpoint=norimagepoint
    imagepoint=image__rationalFunction(hpair,norpoint,nullitytolerance)
    norimagepoint=homogeneousNormalization_pairofcomplexexNullity(imagepoint,
            nullitytolerance)
    number=number+1
    end
    return norimagepoint
end

image__rationalFunction(hpaircompila,(0.5+0.2im , complex(1.0)),15)

abs_image__rationalFunction(hpaircompila,(0.5+0.2im , complex(1.0)),15)

iterationFinalPoint__hpoint(hpaircompila,17,2222,(0.5+0.2im , complex(1.0)))



function valuephimiteration_tolerances_m_hpoint(hpair,
                    nullitytolerance::Int64,
                    cauchytolerance::Int64,
                    seconditernumber::Int64,
                    hpoint::Tuple{Complex{Float64},Complex{Float64}})
    tolnull=1.0/10^(nullitytolerance)
    tolcauchy=1.0/10^(cauchytolerance)
    #zerozero=(complex(0.0),complex(0.0))
    number = 0
    abshpoint=(abs(hpoint[1])+abs(hpoint[2]))
    norpoint=homogeneousNormalization_pairofcomplexexNullity(hpoint,
        nullitytolerance)
    imagepoint=image__rationalFunction(hpair,norpoint,nullitytolerance)
    absimagepoint=(abs(imagepoint[1])+abs(imagepoint[2]))
    norimagepoint=homogeneousNormalization_pairofcomplexexNullity(imagepoint,
        nullitytolerance)
    diffabs=abs(absimagepoint-abshpoint)
    while (number < seconditernumber) & (diffabs>tolcauchy)
    abshpoint=(abs(imagepoint[1])+abs(imagepoint[2]))
    norpoint=norimagepoint
    imagepoint=image__rationalFunction(hpair,norpoint,nullitytolerance)
    norimagepoint=homogeneousNormalization_pairofcomplexexNullity(imagepoint,
            nullitytolerance)
    absimagepoint=(abs(imagepoint[1])+abs(imagepoint[2]))
    diffabs=abs(absimagepoint-abshpoint)
    number=number+1
    end
    return absimagepoint, number
end

valuephimiteration_tolerances_m_hpoint(hpaircompila,45,13,540,
    (0.5+0.2im , complex(1.0)))

function limitphiorbit_iterationvalue__hpair_tolerances_nm_hpoint(hpair,
     nullitytolerance::Int64,
     cauchytolerance::Int64,
     firtsiternumber::Int64,
     seconditernumber::Int64,
     hpoint::Tuple{Complex{Float64},Complex{Float64}})
     hpoint1=iterationFinalPoint__hpoint(hpair, nullitytolerance,
            firtsiternumber,hpoint)
     res=valuephimiteration_tolerances_m_hpoint(hpair,
        nullitytolerance,cauchytolerance,
        seconditernumber,hpoint1)
     return res
end

realinteger1=limitphiorbit_iterationvalue__hpair_tolerances_nm_hpoint(
    hpaircompila,15,25,500,1500, (0.5+0.2im , complex(1.0)))

function threepositioniterations(realinteger, flotanttolerace, maxnumber)
    if realinteger[2]==maxnumber
        pos=0
    elseif realinteger[1]<flotanttolerace
        pos=1
    else
        pos=2
    end
    return pos, realinteger[2]
end
#threepositioniterations(par1 ,   1.0/10^15, 500)

function twoarray_threepositioniterations(hpair,flotanttolerace,
    nullitytolerance::Int64,cauchytolerance::Int64,
    firtsiternumber::Int64, seconditernumber::Int64, rectangle)
    ab=size(rectangle)
    twoarray=[threepositioniterations(limitphiorbit_iterationvalue__hpair_tolerances_nm_hpoint(hpair,
            nullitytolerance,
            cauchytolerance,
            firtsiternumber,
            seconditernumber,
            rectangle[i,j]), flotanttolerace, seconditernumber) for i in 1:ab[1], j in 1:ab[2]]
    return twoarray
end



# cauchytolerance1=15
# flotanttolerace1=1.0/10^cauchytolerance1
# maxnumber1=100
threepositioniterations(realinteger1, 1/10^15,100)


twoarray_threepositioniterations(hpaircompila,10^(-15), 10, 13, 1, 100, rectanglecompila)




function rgb2array_of2arraypositer(colorstrategyname, auxread,
            posread, iterread;
            Nseed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
    plus(threetuple1, threetuple2)=
    (threetuple1[1]+threetuple2[1],threetuple1[2]+threetuple2[2],
        threetuple1[3]+threetuple2[3])
    timestuple(threetuple1, threetuple2)=
    (threetuple1[1]*threetuple2[1],threetuple1[2]*threetuple2[2],
        threetuple1[3]*threetuple2[3])
    leftaction(alpha, threetuple)=
    (alpha*threetuple[1],alpha*threetuple[2],alpha*threetuple[3])
    if colorstrategyname=="position"
    posseed=constructseed(convert(Bool,auxread[1]),convert(Bool,auxread[2]),
            auxread[3];
    Nseed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
    ab=size(posread)
    twoarrayofthreetuples=[posseed[1][posread[i,j]+1]
            for i in 1:ab[1],j in 1:ab[2]]
    elseif  colorstrategyname=="iteration"
    iterseed=constructseed(convert(Bool,auxread[1]),convert(Bool,auxread[2]),
            auxread[4];
    Nseed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
    ab=size(iterread)
    twoarrayofthreetuples=[iterseed[1][iterread[i,j]+1]
            for i in 1:ab[1],j in 1:ab[2]]
    elseif  colorstrategyname=="positioniteration"
    posseed=constructseed(convert(Bool,auxread[1]),convert(Bool,auxread[2]),
            auxread[3];
    Nseed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
    ab=size(posread)
    twoarrayofthreetuples=[
    leftaction(
    (1.0-(iterread[i,j]/(auxread[4]))),
    posseed[1][posread[i,j]+1])
    for i in 1:ab[1], j in 1:ab[2]]
    elseif  colorstrategyname=="iterationposition"
    iterseed=constructseed(convert(Bool,auxread[1]),convert(Bool,auxread[2]),
            auxread[4];
    Nseed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
    iterpos=constructseed(convert(Bool,auxread[1]),convert(Bool,auxread[2]),
            auxread[4];
    Nseed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
    ab=size(iterread)
    len=length(iterpos)
    twoarrayofthreetuples=[
    leftaction(
    (1.0-(posread[i,j]/len)),
    iterseed[1][iterread[i,j]+1])
    for i in 1:ab[1], j in 1:ab[2]]
    elseif  colorstrategyname=="positeraesthetic1"
    posseed=constructseed(convert(Bool,auxread[1]),convert(Bool,auxread[2]),
            auxread[3]+1;
    Nseed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
    ab=size(posread)
    twoarrayofthreetuples=[plus(
    leftaction(
    (1.0-(iterread[i,j]/(auxread[4]))),
    posseed[1][posread[i,j]+1]),
    leftaction((iterread[i,j]/(auxread[4])), posseed[1][posread[i,j]+2]))
    for i in 1:ab[1], j in 1:ab[2]]
    elseif  colorstrategyname=="iterposaesthetic1"
    iterseed=constructseed(convert(Bool,auxread[1]),convert(Bool,auxread[2]),
            auxread[4]+1;
    Nseed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
    ab=size(iterread)
    iterpos=constructseed(convert(Bool,auxread[1]),convert(Bool,auxread[2]),
            auxread[3];
    Nseed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
    len=length(iterpos[1])
    twoarrayofthreetuples=[plus(
    leftaction((1.0-(posread[i,j]/len)),
    iterseed[1][iterread[i,j]+1]),
    leftaction(posread[i,j]/len, iterseed[1][iterread[i,j]+2]))
    for i in 1:ab[1], j in 1:ab[2]]
    end
    s=size(twoarrayofthreetuples)[1]
    twoarrayofthreetuplesnewjulia11=[twoarrayofthreetuples[i,j][k]
        for i in 1:s,  j in 1:s, k in 1:3]
    return  twoarrayofthreetuplesnewjulia11
end

function constructseed(issuperzero::Bool,
                    isinfinity::Bool,
                    numberofordinaryfixedpoint::Int64;
        Nseed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
        if  issuperzero==false && isinfinity==false
        seed1=[Colors.RGB(0.0,0.0,0.0)];numberofcolors=numberofordinaryfixedpoint+1
        elseif  issuperzero==false && isinfinity==true
        seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(1.0,1.0,0.0)]
        numberofcolors=numberofordinaryfixedpoint+2
        elseif  issuperzero==true && isinfinity==false
        seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5)]
        numberofcolors=numberofordinaryfixedpoint+2
        else
        seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5),Colors.RGB(1.0,1.0,0.0)]
        numberofcolors=numberofordinaryfixedpoint+3
        end
    seed=union(seed1,Nseed)
    gseed=Colors.distinguishable_colors(numberofcolors, seed)
    gseedtriple=[(gseed[i].r,gseed[i].g,gseed[i].b) for i in 1:length(gseed)]
    return gseedtriple , gseed
end

function pathsconstructionrgbimages(colorstrategyname;prefix="", extension="")
            if colorstrategyname=="position"
            paths=(joinpath("../images/",  join([prefix,"pos",extension])))
            elseif colorstrategyname=="iteration"
            paths=(joinpath("../images/",  join([prefix,"iter",extension])))
            elseif colorstrategyname=="positioniteration"
            paths=(joinpath("../images/",
                join([prefix,"positer",extension])))
            elseif colorstrategyname=="iterationposition"
            paths=(joinpath("../images/",
                join([prefix,"iterpos",extension])))
            elseif colorstrategyname=="positeraesthetic1"
            paths=(joinpath("../images/",
                join([prefix,"positeraesth1",extension])))
            elseif colorstrategyname=="iterposaesthetic1"
            paths=(joinpath("../images/",
                join([prefix,"iterposaesth1",extension])))
            end
            return paths
end

pathsconstructionrgbimages("position";prefix="", extension=".png")





function rgb2arrayBasinsUnivariateRationalFunctions(
        coefficientlistnum::Array{T,1},
        coefficientlistden::Array{T,1},
        expresolution::Int=8,
        iterationmax::Int=25,
        tolerances::Tuple{Int64,Int64,Int64}=(15,3,3);
        colorstrategy::AbstractString="positioniteration",
        model::AbstractString="localrectangle",
        rectanglesides::Tuple{Float64,Float64,Float64,Float64}=
            (-1.5,1.5,-1.5,1.5)) where {T<:Number}
            modeldom=eval(model)
            rectanglesidesmod=rectanglesides
            colorstrategymod=colorstrategy
        newcoefficientlist=samenumberofcoeficcients(coefficientlistnum,coefficientlistden)
        fixedpointsmatrixpairpositerex=fixedPointListex_matrixpositioninterations_RationalFunction(newcoefficientlist[1],
        newcoefficientlist[2], expresolution,iterationmax,
        tolerances,
         model,rectanglesides)
         auxfixindmatrixpi=auxfixedPointListex_matrixpositioninterations_RationalFunction(
                 coefficientlistnum,
                 coefficientlistden,
                 expresolution,
                 iterationmax,
                 tolerances,
                 modeldom,
                 rectanglesidesmod)
                 ab=size(auxfixindmatrixpi[2])
                  posred=[auxfixindmatrixpi[2][i,j][1] for i in 1:ab[1] , j in 1:ab[2]]
                  iterread=[auxfixindmatrixpi[2][i,j][2] for i in 1:ab[1] , j in 1:ab[2]]
                  rectanglesidesmod=[rectanglesidesmod[1],rectanglesidesmod[2],rectanglesidesmod[3],rectanglesidesmod[4]]
                 rgb2array1=rgb2array_of2arraypositer(colorstrategymod, auxfixindmatrixpi[1],
                            posred, iterread,
                             Nseed=[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
                #plotresult=PyPlot.imshow(twoarray3tuplergb;extent=rectanglesidesmod)
        return   rgb2array1
end


rgb2arrayBasinsUnivariateRationalFunctions(rfcompila[1],rfcompila[2],3)



function rgb2arrayBasinsFixIndUnivariateRationalFunctions(
        coefficientlistnum::Array{T,1},
        coefficientlistden::Array{T,1},
        expresolution::Int=8,
        iterationinterval::Tuple{Int,Int}=(1,100),
        tolerances::Tuple{Int64,Int64,Int64}=(15,3,3);
        colorstrategy::AbstractString="positioniteration",
        model::AbstractString="localrectangle",
        rectanglesides::Vector{Float64}=
            (-1.5,1.5,-1.5,1.5)) where {T<:Number}
            modeldomain=model
            rectanglesidesmod=rectanglesides
            colorstrategymod=colorstrategy
        newcoefficientlist=samenumberofcoeficcients(coefficientlistnum,coefficientlistden)
        hpair1=pairoffunctions(newcoefficientlist[1], newcoefficientlist[2])
        nullitytolerance1=tolerances[1]
        cauchytolerance1=tolerances[2]
        flotanttolerace1=1.0/10^cauchytolerance1
        firtsiternumber1=iterationinterval[1]
        seconditernumber1=iterationinterval[2]
        xinterv=(rectanglesides[1],rectanglesides[2])
        yinterv=(rectanglesides[3],rectanglesides[4])
        if modeldomain=="localrectangle"
            rect=rectangle(xinterv,  yinterv, expresolution)
        elseif modeldomain=="invertedlocalrectangle"
            rect=invertedrectangle(xinterv,  yinterv, expresolution)
        end
        tatp=twoarray_threepositioniterations(hpair1,
        flotanttolerace1,nullitytolerance1,cauchytolerance1,firtsiternumber1,
        seconditernumber1, rect)
        rcg=[Colors.RGB(1.,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)]
        auxnread=(0,0,3,seconditernumber1)
        ab=size(tatp)
        posnread=[tatp[i,j][1] for i in 1:ab[1], j in 1:ab[2]]
        iternread=[tatp[i,j][2] for i in 1:ab[1], j in 1:ab[2]]
        prergb=rgb2array_of2arraypositer(colorstrategy, auxnread,
            posnread, iternread;Nseed=rcg)
        return  prergb
end




function xyzlimitphiorbit(coefficientlistnum,coefficientlistden,
    nullitytolerance::Int64=10,
    cauchytolerance::Int64=15,
    firtsiternumber::Int64=100,
    seconditernumber::Int64=100,
    expresolution::Int64=8,rectanglevector::Vector{Float64}=[-1.5,1.5,-1.5,1-5])
    hpair1=pairoffunctions(coefficientlistnum,coefficientlistden)#ho=println("hola")
    numberres=2^expresolution
    x=range(rectanglevector[1], stop = rectanglevector[2], length = numberres)
    y=range(rectanglevector[3], stop = rectanglevector[4], length = numberres)
    z=[limitphiorbit_iterationvalue__hpair_tolerances_nm_hpoint(
            hpair1,
            nullitytolerance,
            cauchytolerance,
            firtsiternumber,
            seconditernumber,
     homogeneousNormalization_pairofcomplexexNullity(
                (complex(x[i]+im*y[j]),complex(1.0)),20))[1] for i=1:length(x) , j=1:length(y)]
    return [x, y, z]
end



function polynomialToRelaxedNewtonRationalFunction(coefficienlist, h)
    pol=Polynomials.Poly(coefficienlist)
    polx=Polynomials.Poly([0.0,1.0])
    numerator=Polynomials.polyder(pol)*polx-h*pol
    denominator=Polynomials.polyder(pol)
    coeffnumerator=complex([numerator[i] for i in 0:1:length(numerator)])
    coeffdenominator=complex([denominator[i] for i in 0:1:length(denominator)])
    nundem=samenumberofcoeficcients(coeffnumerator, coeffdenominator)
    nundemnorm=normalize_pairhpolynomials(nundem[1],nundem[2])
    return  nundemnorm
end


end # module
