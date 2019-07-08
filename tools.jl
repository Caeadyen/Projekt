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
  end
  for i = 1:2:(2*parameter.N1+1)
    plot(data.x_node[i:2*parameter.N1+1:end,1],data.x_node[i:2*parameter.N1+1:end,2],"b-")
  end
end

function plotpressure(data,parameter)
  x_node=Int64[]

  for i = (2*parameter.N2+1):-2:1
      for j=1+((i-1)*(2*parameter.N1+1)):2:(i*(2*parameter.N1+1))
      push!(x_node,j)
    end
  end

  X = reshape(data.x_node[x_node,1], parameter.N1+1,parameter.N2+1)'
  Y = reshape(data.x_node[x_node,2], parameter.N1+1,parameter.N2+1)'
  p=reshape(data.p',parameter.N1+1,parameter.N2+1)'* 1e-9
  contourf(X,Y,p)
  colorbar()

end
