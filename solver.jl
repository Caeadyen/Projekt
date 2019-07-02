#assemble Steifigkeitsmatrix und Lastvektor fuer ein finites Element
function assemble_element(i,data,parameter)
  quad = data.square[i,:]
  corner = data.corner[i,:]
  points, weights = gausslegendrerule(3)
  quad_pairs = zip(points,weights)
  Kvve=zeros(length(quad),length(quad))
  Kvpe=zeros(length(quad),length(corner))
  Kppe=zeros(length(corner),length(corner))
  Fk =zeros(1,length(quad))
  #TODO nu,E in parameter und clean up was wo berechnet wird
  #TODO J detJ J^-1 einbaun
  #D = [2 0 0; 0 2 0; 0 0 1]
  nu=0.35
  E=3*10^10
  Ds=(E/(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2]
  D=[4/3 -2/3 0;-2/3 4/3 0;0 0 1]
  G=E/(2*(1+nu))
  dt=0.1*parameter.my/G
println(dt)
  dN=zeros(2,length(quad))
  dNdX=zeros(2,length(quad))
  B=zeros(3,length(quad))
  Bvol=zeros(1,length(quad))
  Fv=zeros(1,length(quad))
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
    dN2,NN2 = shape(x1,x2,1)
    jac=data.x_node[quad[1:2:end],:]'*dN'
    detJ=jac[1,1]*jac[2,2]-jac[1,2]*jac[2,1]
    invj=inv(jac)
    #println(detJ)
    dNdX=dN'*invj
    dNdX=dNdX'
    mu=1/((1/parameter.my)+(1/(G*dt)))
    chi = 1/(1+G*dt/parameter.my)
    B[2,indexe]=dNdX[2,:]
    B[3,indexe]=dNdX[1,:]
    B[1,indexu]=dNdX[1,:]
    B[3,indexu]=dNdX[2,:]
    tau=Ds*B*data.u[quad]
    Bvol[1,indexu]=dNdX[1,:]
    Bvol[1,indexe]=dNdX[2,:]
    #W= [0 0 2; 0 0 -2; -1 1 0]*Bw*data.u[quad]*dt

    Kvve += w1*w2*detJ*mu*((B'*D)*B)
    Kvpe += -w1*w2*(Bvol'*NN2)*detJ

    Kppe += w1*w2*detJ*(NN2'*NN2)
    x=w1*w2*detJ*dNdX'*[tau[1] ; tau[3]]*chi
    y=w1*w2*detJ*dNdX'*[tau[3] ; tau[2]]*chi
    Fv[indexu]=x
    Fv[indexe]=y
    Fk += ((w1*w2*detJ)*parameter.f*NN) -Fv + (parameter.boundary*data.n_boundary[quad,:]' .*(w1*w2*detJ*[1 0]*NN))
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
    #Sk=zeros(18,18)
    Fk =zeros(1,18)
    Kvve, Kvpe ,Kppe, Fk = assemble_element(i,data,parameter)
    quad = data.square[i,:]
    corner = data.corner[i,:]
    F[quad]=F[quad]+transpose(Fk)

    for j = 1:18
            for i = 1:18
                    push!(Ivv, quad[i])
                    push!(Jvv, quad[j])
                    push!(Vvv, Kvve[i,j])

            end
    end


    for j = 1:18
            for i = 1:4
                    push!(Ivp, corner[i])
                    push!(Jvp, quad[j])
                    push!(Vvp, Kvpe[j,i])
            end
    end


  for j = 1:4
          for i = 1:4
                  push!(Ipp, corner[i])
                  push!(Jpp, corner[j])
                  push!(Vpp, Kppe[j,i])
          end
  end
end

  Kvv=sparse(Ivv,Jvv,Vvv,data.n_nd,data.n_nd)
  Kvp=sparse(Ivp,Jvp,Vvp,data.n_corner,data.n_nd)
  Kpp=sparse(Ipp,Jpp,Vpp,data.n_corner,data.n_corner)

  return Kvv,Kvp,Kpp,F

end
