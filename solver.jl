#assemble Steifigkeitsmatrix und Lastvektor fuer ein finites Element
function assemble_element(i,data,parameter)
  quad = data.square[i,:]
  corner = data.corner[i,:]
  points, weights = gausslegendrerule(3)
  quad_pairs = zip(points,weights)
  Kvve=zeros(length(quad),length(quad))
  Kvpe=zeros(length(quad),length(corner))
  Fk =zeros(1,length(quad))
  #TODO nu,E in parameter und clean up was wo berechnet wird
  D = [2 0 0; 0 2 0; 0 0 1]
  nu=0.35
  E=3*10^10
  EM=(E/(1-nu^2))*[1 nu 0;nu 1 0;0 0 (1-nu)/2]
  G=E/(2*(1+nu))
  dt=parameter.my/G
  dN=zeros(2,length(quad))
  B=zeros(3,length(quad))
  Bw=zeros(2,length(quad))
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
    B[2,indexe]=dN[2,:]
    B[3,indexe]=dN[1,:]
    B[1,indexu]=dN[1,:]
    B[3,indexu]=dN[2,:]
    tau=E*B*data.u[quad]

    #W= [0 0 2; 0 0 -2; -1 1 0]*Bw*data.u[quad]*dt
    Kvve += w1*w2*J*((B'*(2*parameter.my*dt/(dt*parameter.my/G))*D)*B)
    Kvpe += -w1*w2*(B[1:2,:]'*[NN2;NN2])*J
    Fk += ((w1*w2*J)*parameter.f*NN) - tau'*B*((parameter.my/G)/(dt*(parameter.my/G)))+ (parameter.boundary*data.n_boundary[quad,:]' .*(w1*w2*J*[1 0]*NN))
    end
return Kvve, Kvpe, Fk
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
    Sk=zeros(18,18)
    Fk =zeros(1,18)
    Kvve, Kvpe ,Fk = assemble_element(i,data,parameter)
    quad = data.square[i,:]
    F[quad]=F[quad]+transpose(Fk)

    for j = 1:18
            for i = 1:18
                    push!(Ivv, quad[i])
                    push!(Jvv, quad[j])
                    push!(Vvv, Kvve[i,j])

            end
    end

    quad = data.square[i,:]
    corner = data.corner[i,:]
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
