#Data structure to store the mesh
mutable struct Data
  dim  :: Int

  n_el :: Int
  n_nd :: Int
  n_corner :: Int


  x_node   :: Array{Float64, 2}
  square   :: Array{Int, 2}

  n_boundary :: Array{Int, 2}

  x_corner   :: Array{Float64, 2}

  corner   :: Array{Int, 2}

  c_boundary :: Array{Int, 2}

  u :: Array{Float64}

  p :: Array{Float64}

  tau :: Array{Float64, 3}

end

#Data structure for in Parameter
#L1/L2 L�nge des gebiets
#N1/N2 in wieviele unterteilt wird
#f (spaltenvektor), Volume force
#sub Anzahl in die das Gebiet in y-Richtung unterteilt werden soll, N2 muss durch diese Zahl teilbar sein

#my (spaltenvektor), Viskosität
#nu Poisson's ratio
#E Young's modulus
#rbounadry Speed the right bounadry moves in m/year
#lbounadry Speed the right bounadry moves in m/year
#dt in years
#maxtime how many years it should run


struct Parameter

    L1:: Int
    L2:: Int
    N1:: Int
    N2:: Int
    f :: Array{Int, 1}
    sub :: Int
    my :: Array{Float64, 1}
    nu :: Float16
    E :: Float64
    G :: Float16
    lboundary :: Float16
    rboundary :: Float16
    dt :: Int
    maxtime :: Int

end



function gausslegendrerule(n::Integer)
  lambda, V = eigen(.5 * SymTridiagonal(ones(n), 1 ./sqrt.(4 .- 1 ./ (1:(n-1)).^2)));
  return lambda, V[1,:].^2
end

function plotgrid(data,parameter)
  for i = 1:2:(2*parameter.N2+1)
    plot(data.x_node[(1+((i-1)*(2*parameter.N1+1))):(i*(2*parameter.N1+1)),1],data.x_node[(1+((i-1)*(2*parameter.N1+1))):(i*(2*parameter.N1+1)),2],"b-")
    # println(data.x_node[(1+((i-1)*(2*parameter.N1+1))):(i*(2*parameter.N1+1)),1])
    # println(data.x_node[(1+((i-1)*(2*parameter.N1+1))):(i*(2*parameter.N1+1)),2])
  end
  for i = 1:2:(2*parameter.N1+1)
    plot(data.x_node[i:2*parameter.N1+1:end,1],data.x_node[i:2*parameter.N1+1:end,2],"b-")
  end
end

function plotpressure(data,parameter)
  x_node=Int64[]

  for i = 1:2:(2*parameter.N2+1)
      for j=1+((i-1)*(2*parameter.N1+1)):2:(i*(2*parameter.N1+1))
      push!(x_node,j)
    end
  end

  X = reshape(data.x_node[x_node,1], parameter.N1+1,parameter.N2+1)'
  Y = reshape(data.x_node[x_node,2], parameter.N1+1,parameter.N2+1)'
  p=reshape(data.p',parameter.N1+1,parameter.N2+1)'*1e-9


  contourf(X,Y,p)
  colorbar()

end

function plottau(data,parameter,tau,option)
  x_node=Int64[]

  for i = 1:2:(2*parameter.N2+1)
      for j=1+((i-1)*(2*parameter.N1+1)):2:(i*(2*parameter.N1+1))
      push!(x_node,j)
    end
  end
  if(option=="tau2nd")
      sigm=1e-9*(0.5*tau[1].^2+0.5*tau[2].^2+tau[3]).^0.5
      X = reshape(data.x_node[x_node,1], parameter.N1+1,parameter.N2+1)'
      Y = reshape(data.x_node[x_node,2], parameter.N1+1,parameter.N2+1)'
      sigm=reshape(sigm,parameter.N1+1,parameter.N2+1)'
      contourf(X,Y,sigm, norm=matplotlib.colors.LogNorm())
      colorbar()
    elseif(option=="tauxx")
        sigm = tau[1]
        X = reshape(data.x_node[x_node,1], parameter.N1+1,parameter.N2+1)'
        Y = reshape(data.x_node[x_node,2], parameter.N1+1,parameter.N2+1)'
        sigm=reshape(sigm,parameter.N1+1,parameter.N2+1)'
        contourf(X,Y,sigm)
        colorbar()
    elseif(option=="tauyy")
        sigm = tau[2]
        X = reshape(data.x_node[x_node,1], parameter.N1+1,parameter.N2+1)'
        Y = reshape(data.x_node[x_node,2], parameter.N1+1,parameter.N2+1)'
        sigm=reshape(sigm,parameter.N1+1,parameter.N2+1)'
        contourf(X,Y,sigm)
        colorbar()
    elseif(option=="tauxy")
        sigm = tau[3]
        X = reshape(data.x_node[x_node,1], parameter.N1+1,parameter.N2+1)'
        Y = reshape(data.x_node[x_node,2], parameter.N1+1,parameter.N2+1)'
        sigm=reshape(sigm,parameter.N1+1,parameter.N2+1)'
        contourf(X,Y,sigm)
        colorbar()
    else
        println("Kannt ", option, " nicht printen")
    end


end

function  convert_tau(data,parameter)
    tau_xx = zeros(1, (parameter.N1+1)*(parameter.N2+1))
    tau_xy = zeros(1, (parameter.N1+1)*(parameter.N2+1))
    tau_yy = zeros(1, (parameter.N1+1)*(parameter.N2+1))
    for j=1:parameter.N2
        for i=1:parameter.N1
        # punkt 1
        tau_xx[(parameter.N1+1)*(j-1)+i] += data.tau[1,1,parameter.N1*(j-1)+i]
        tau_xy[(parameter.N1+1)*(j-1)+i] += data.tau[3,1,parameter.N1*(j-1)+i]
        tau_yy[(parameter.N1+1)*(j-1)+i] += data.tau[2,1,parameter.N1*(j-1)+i]
        # punkt 3
        tau_xx[(parameter.N1+1)*(j-1)+i+parameter.N1+1] += data.tau[1,2,parameter.N1*(j-1)+i]
        tau_xy[(parameter.N1+1)*(j-1)+i+parameter.N1+1] += data.tau[3,2,parameter.N1*(j-1)+i]
        tau_yy[(parameter.N1+1)*(j-1)+i+parameter.N1+1] += data.tau[2,2,parameter.N1*(j-1)+i]
        # punkt 7
        tau_xx[(parameter.N1+1)*(j-1)+i+1] += data.tau[1,3,parameter.N1*(j-1)+i]
        tau_xy[(parameter.N1+1)*(j-1)+i+1] += data.tau[3,3,parameter.N1*(j-1)+i]
        tau_yy[(parameter.N1+1)*(j-1)+i+1] += data.tau[2,3,parameter.N1*(j-1)+i]
        # punkt 9
        tau_xx[(parameter.N1+1)*(j-1)+i+parameter.N1+2] += data.tau[1,4,parameter.N1*(j-1)+i]
        tau_xy[(parameter.N1+1)*(j-1)+i+parameter.N1+2] += data.tau[3,4,parameter.N1*(j-1)+i]
        tau_yy[(parameter.N1+1)*(j-1)+i+parameter.N1+2] += data.tau[2,4,parameter.N1*(j-1)+i]
    end
end

    debug_c1 = 0
    debug_c2 = 0
    debug_c3 = 0
    for i=1:((parameter.N1+1)*(parameter.N2+1))
        if sum(data.c_boundary[i,:])==0
            tau_xx[i] /= 4
            tau_xy[i] /= 4
            tau_yy[i] /= 4
            debug_c1 += 1
        end
        if sum(data.c_boundary[i,:])==1
            tau_xx[i] /= 2
            tau_xy[i] /= 2
            tau_yy[i] /= 2
            debug_c2 += 1
        end
        if sum(data.c_boundary[i,:])==2
            tau_xx[i] /= 1
            tau_xy[i] /= 1
            tau_yy[i] /= 1
            debug_c3 += 1
        end
    end
    return tau_xx,tau_xy,tau_yy
end
