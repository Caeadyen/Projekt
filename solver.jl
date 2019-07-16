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

  Ds=(parameter.E/(1-parameter.nu^2))*[1 parameter.nu 0;parameter.nu 1 0;0 0 (1-parameter.nu)/2]
  D=[4/3 -2/3 0;-2/3 4/3 0;0 0 1]

  dt=3600*24*365.25*parameter.dt #100jahre
  dN=zeros(2,length(quad))
  dNdX=zeros(2,length(quad))
  B=zeros(3,length(quad))
  Bw=zeros(1,length(quad))
  Bvol=zeros(1,length(quad))
  N=zeros(2,convert(Int,(length(quad)/2)))
  indexu=Int64[]
  indexe=Int64[]



  for i = 1:convert(Int,(length(quad)/2))
    push!(indexu,2*i-1)
    push!(indexe,2*i)
  end

  count = 1
  count_tau = 1

  for (x1,w1) in quad_pairs, (x2,w2) in quad_pairs

    dN,NN = shape(x1,x2,2)
    dN2,NN2 = shape(x1,x2,1)
    #Jacobian, its Determinant and inverse
    jac=data.x_node[quad[1:9],:]'*dN'
    detJ=jac[1,1]*jac[2,2]-jac[1,2]*jac[2,1]
    invj=(1/detJ).*[jac[2,2] -jac[1,2]; -jac[2,1] jac[1,1]]

    #derivatives
    dNdX=dN'*invj

    #B,W Matrices
    B[2,10:18]=dNdX[:,2]
    B[3,1:9]=dNdX[:,2]
    B[1,1:9]=dNdX[:,1]
    B[3,10:18]=dNdX[:,1]
    Bw[1,1:9]=-0.5*dNdX[:,2]
    Bw[1,10:18]=0.5*dNdX[:,1]
    Bvol[1,1:9]=dNdX[:,1]
    Bvol[1,10:18]=dNdX[:,2]
    W=[0 0 2; 0 0 -2;-1 1 0].*(Bw*data.u[quad])



    for j = 1:parameter.sub
       # Jeder Schicht wird das richtig my zugewiesen, d.h. den unteren N1 Elementen das erste my, den darüberliegenden N1 Elemente das zweite my usw.
        if(((i-1)÷(convert(Int,(data.n_el/parameter.sub)))+1) == j)
            mu=1/((1/parameter.my[j])+(1/(parameter.G*dt)))
            chi =1/(1+(parameter.G*dt/parameter.my[j]))
            #Stress vector
            tau=mu*Ds*B*(data.u[quad])
            #saving stress near the 4 corners for plotting
            if(count == 1 || count == 3 || count == 7 || count == 9)
              data.tau[:,count_tau, i] = tau
              count_tau += 1
            end

            Kvve += w1*w2*detJ*mu*((B'*D)*B)
            Kvpe += -w1*w2*(Bvol'*NN2)*detJ
            Fk += [ (w1*w2*detJ)*parameter.f[1]*(NN[1,:].*NN[2,:]);(w1*w2*detJ)*parameter.f[2]*(NN[1,:].*NN[2,:])] -w1*w2*detJ*B'*(I+W)*tau*chi
            count += 1
        end
    end
  end
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
