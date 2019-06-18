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
  #Mesh(dim::Int, n_el::Int, n_nd::Int, n_corner :: Int, x_node::Array{Float64,2}, square::Array{Int64,2}, n_boundary::Array{Int64,2}, x_corner::Array{Float64, 2},corner::Array{Int64,2}, c_boundary::Array{Int64,2}) = new(dim, n_el, n_nd, n_corner,x_node, square, n_boundary, x_corner, corner, c_boundary)
end

#Data structure for in Parameter
#L1/L2 Länge des gebiets
#N1/N2 in wieviele unterteilt wird
#Boundary grenzen [linke Seite, untere Seite, rechte Seite obere Seite]
#f
#my
#k
struct Parameter

    L1:: Int
    L2:: Int
    N1:: Int
    N2:: Int
    boundary :: Array{Int, 2}
    f :: Array{Int, 2}
    my :: Int
    k :: Float16

end



function gausslegendrerule(n::Integer)
  λ, V = eigen(.5 * SymTridiagonal(ones(n), 1 ./sqrt.(4 .- 1 ./ (1:(n-1)).^2)));
  return λ, V[1,:].^2
end

function plotgrid(data,parameter)
  for i = 1:parameter.N2+1
    plot(data.x_corner[(1+((i-1)*(parameter.N1+1))):(i*(parameter.N1+1)),1],data.x_corner[(1+((i-1)*(parameter.N1+1))):(i*(parameter.N1+1)),2],"b-")
  end
  for i = 1:parameter.N1+1
    plot(data.x_corner[i:(parameter.N1+1):end,1],data.x_corner[i:(parameter.N1+1):end,2],"b-")
  end
end
