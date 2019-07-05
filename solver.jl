#assemble Steifigkeitsmatrix und Lastvektor fuer ein finites Element
function assemble_element(i,data,parameter)
  quad = data.square[i,:]
  corner = data.corner[i,:]
  points, weights = gausslegendrerule(3)
  quad_pairs = zip(points,weights)
  Kvve=zeros(length(quad),length(quad))
  Kvpe=zeros(length(quad),length(corner))
  Kppe=zeros(length(corner),length(corner))
  Fk =zeros(length(quad))
  Ft =zeros(2,9)

  nu=parameter.nu
  E=parameter.E
  Ds=(E/(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2]
  D=[4/3 -2/3 0;-2/3 4/3 0;0 0 1]

  G=1e11
  dt=3600*24*365.25 #100jahre
  dN=zeros(2,length(quad))
  dNdX=zeros(2,length(quad))
  B=zeros(3,length(quad))
  Bvol=zeros(1,length(quad))
  Fv=zeros(1,length(quad))
  N=zeros(2,convert(Int,(length(quad)/2)))
  indexu=Int64[]
  indexe=Int64[]


  for i = 1:convert(Int,(length(quad)/2))
    push!(indexu,2*i-1)
    push!(indexe,2*i)
  end


  for (x1,w1) in quad_pairs, (x2,w2) in quad_pairs
    dN,NN = shape(x1,x2,2)
    dN2,NN2 = shape(x1,x2,1)
    jac=data.x_node[quad[1:9],:]'*dN'
    detJ=jac[1,1]*jac[2,2]-jac[1,2]*jac[2,1]
    invj=inv(jac)
    dNdX=dN'*invj
    mu=1/((1/parameter.my)+(1/(G*dt)))
    chi =1/(1+(G*dt/parameter.my))
    B[2,10:18]=dNdX[:,2]
    B[3,1:9]=dNdX[:,2]
    B[1,1:9]=dNdX[:,1]
    B[3,10:18]=dNdX[:,1]

    tau=Ds*B*(data.u[quad])
    t=[tau[1] tau[3];tau[3] tau[2]]
    Bvol[1,1:9]=dNdX[:,1]
    Bvol[1,10:18]=dNdX[:,2]

    #W= [0 0 2; 0 0 -2; -1 1 0]*Bw*data.u[quad]*dt

    Kvve += w1*w2*detJ*mu*((B'*D)*B)
    Kvpe += -w1*w2*(Bvol'*NN2)*detJ
    Fv=w1*w2*detJ*t*dNdX'*chi
    Ft += ((w1*w2*detJ)*parameter.f*(NN[1,:]'.*NN[2,:]')) -Fv

  end

  Fk=[Ft[1,:] ; Ft[2,:]]


return Kvve, Kvpe,Fk
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

  for i = 1:data.n_el

    Fk =zeros(1,18)
    Kvve, Kvpe , Fk = assemble_element(i,data,parameter)
    quad = data.square[i,:]
    corner = data.corner[i,:]
    F[quad]=F[quad]+Fk
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

end

  Kvv=sparse(Ivv,Jvv,Vvv,data.n_nd,data.n_nd)
  Kvp=sparse(Ivp,Jvp,Vvp,data.n_corner,data.n_nd)

  return Kvv,Kvp,F

end
