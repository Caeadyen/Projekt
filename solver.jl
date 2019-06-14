#assemble Steifigkeitsmatrix und Lastvektor fuer ein finites Element
function assemble_element(i,data,parameter)
  quad = data.square[i,:]
  points, weights = gausslegendrerule(3)
  quad_pairs = zip(points,weights)
  Kvve=zeros(length(quad),length(quad))
  Fk =zeros(1,length(quad))
  D = [4/3 -2/3 0; -2/3 4/3 0; 0 0 1]
  dN=zeros(2,length(quad))
  B=zeros(3,length(quad))
  N=zeros(2,convert(Int,(length(quad)/2)))
  indexu=Int64[]
  indexe=Int64[]
  J= 0.5*abs((data.x_node[quad[1],2]-data.x_node[quad[5],2])*(data.x_node[quad[7],1]-data.x_node[quad[3],1])+(data.x_node[quad[3],2]-data.x_node[quad[7],2])*(data.x_node[quad[1],1]-data.x_node[quad[5],1]))

  for i = 1:convert(Int,(length(quad)/2))
    push!(indexu,2*i-1)
    push!(indexe,2*i)
  end

  for (x1,w1) in quad_pairs, (x2,w2) in quad_pairs
    dN,NN = shape(x1,x2,2)
    B[2,indexe]=dN[2,:]
    B[3,indexe]=dN[1,:]
    B[1,indexu]=dN[1,:]
    B[3,indexu]=dN[2,:]

    Kvve += w1*w2*J*((B'*D)*B)*parameter.my
    #F berechnen und für links und rechts das drücken draufrechnen
    #TODO f,t in parameter struct
    Fk += ((w1*w2*J)*[0 0]*NN) + (data.n_boundary[quad,:]*parameter.boundary' .*((w1*w2*J*[1 0])*NN)')'

  end

  corner = data.corner[i,:]
  points, weights = gausslegendrerule(2)
  quad_pairs = zip(points,weights)
  Kvpe=zeros(length(quad),length(corner))
  Kppe=zeros(length(corner),length(corner))
  indexu=Int64[]
  indexe=Int64[]
#building des index vectors, ging irgendwie nicht mit collect
  for i = 1:convert(Int,(length(corner)/2))
    push!(indexu,2*i-1)
    push!(indexe,2*i)
  end
  for (x1,w1) in quad_pairs, (x2,w2) in quad_pairs
    dN,NN = shape(x1,x2,1)
    Kvpe += w1*w2*(B[1:2,:]'*dN)*J
    Kppe += w1*w2*dN'*dN
  end

return Kvve, Kvpe, Kppe, Fk
end

#assemble gesamte Steifigkeitsmatrix und Lastvektor
function assemble_whole(data,parameter)

  F = zeros(data.n_nd)

  Ivv=Int64[]
  Jvv=Int64[]
  Vvv=Float64[]

  Ivp=Int64[]
  Jvp=Int64[]
  Vvp=Float64[]

  Ipp=Int64[]
  Jpp=Int64[]
  Vpp=Float64[]

  for i = 1:data.n_el
    Sk=zeros(18,18)
    Fk =zeros(1,18)
    Kvve, Kvpe ,Kppe , Fk = assemble_element(i,data,parameter)
    quad = data.square[i,:]
    F[quad]=F[quad]+transpose(Fk)

    for j = 1:18
            for i = 1:18
                    push!(Ivv, quad[i])
                    push!(Jvv, quad[j])
                    push!(Vvv, Kvve[i,j])

            end
    end

    quad = data.square[i,[1,3,5,7]]
    corner = data.corner[i,:]
    for j = 1:4
            for i = 1:4
                    push!(Ivp, corner[i])
                    push!(Jvp, quad[j])
                    push!(Vvp, Kvpe[i,j])
            end
    end



    corner = data.corner[i,:]
    for j = 1:4
          for i = 1:4
                push!(Ipp, corner[i])
                push!(Jpp, corner[i])
                push!(Vpp, Kppe[i,j])
          end
    end
  end


  Kvv=sparse(Ivv,Jvv,Vvv,data.n_nd,data.n_nd)
  Kvp=sparse(Ivp,Jvp,Vvp,data.n_corner,data.n_nd)
  Kpp=sparse(Ipp,Jpp,Vpp,data.n_corner,data.n_corner)

  return Kvv,Kvp,Kpp,F

end
