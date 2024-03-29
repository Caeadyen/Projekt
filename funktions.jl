using SparseArrays
using LinearAlgebra
include("tools.jl")
#Building the grid
function builddata(Parameter)

  h1=Parameter.L1/(2*Parameter.N1)
  h2=Parameter.L2/(2*Parameter.N2)
  dim = 2
  n_el = Parameter.N1*Parameter.N2
  n_nd = ((2*Parameter.N1)+1)*((2*Parameter.N2)+1)
  n_corner = (Parameter.N1+1)*(Parameter.N2+1)
  n_boundary = zeros(Int,n_nd*2,4)
  x_node = zeros(Float64,n_nd,2)
  c_boundary = zeros(Int,n_corner,4)
  x_corner = zeros(Float64,n_corner,2)
  square = zeros(Int,n_el,18)
  corner = zeros(Int,n_el,4)
  n_count = 1

  for j = 1:((2*Parameter.N2)+1)
    yj = (j-1)*h2
    for i = 1:((2*Parameter.N1)+1)
      xi = (i-1)*h1
      x_node[n_count,1] = xi
      x_node[n_count,2] = yj
      #festlegen der Boundary: linker Rand= 1; unterer =2; rechter = 3; oberer = 4
      if(i == 1)
        n_boundary[n_count,1] = 1
        n_boundary[n_count+n_nd,1] = 1
      end
      if(i==(2*Parameter.N1)+1)
        n_boundary[n_count,3] = 1
        n_boundary[n_count+n_nd,3] = 1
      end
      if(j==1)
          n_boundary[n_count,2] = 1
          n_boundary[n_count+n_nd,2] = 1
      end
      if(j==(2*Parameter.N2)+1)
          n_boundary[n_count,4] = 1
          n_boundary[n_count+n_nd,4] = 1
      end
      n_count += 1
    end
  end

  e_count = 0

  for j = 2:2:2*Parameter.N2
    for i = 1:Parameter.N1
      e_count += 1
      square[e_count,1] = (j-2)* (2*Parameter.N1+1) + (2*i)-1
      square[e_count,10] = square[e_count,1]+n_nd

      square[e_count,2] = (j-2)*(2*Parameter.N1+1) + (2*i)+1
      square[e_count,11] = square[e_count,2]+n_nd

      square[e_count,3] = (j)*(2*Parameter.N1+1) + (2*i)+1
      square[e_count,12] = square[e_count,3]+n_nd

      square[e_count,4] = (j)*(2*Parameter.N1+1) + (2*i)-1
      square[e_count,13] = square[e_count,4]+n_nd

      square[e_count,5] = (j-2)* (2*Parameter.N1+1) + (2*i)
      square[e_count,14] = square[e_count,5]+n_nd

      square[e_count,6] = (j-1)* (2*Parameter.N1+1) + (2*i)+1
      square[e_count,15] = square[e_count,6]+n_nd

      square[e_count,7] = (j)* (2*Parameter.N1+1) + (2*i)
      square[e_count,16] = square[e_count,7]+n_nd

      square[e_count,8] = (j-1)* (2*Parameter.N1+1) + (2*i)-1
      square[e_count,17] = square[e_count,8]+n_nd

      square[e_count,9] = (j-1)* (2*Parameter.N1+1) + (2*i)
      square[e_count,18] = square[e_count,9]+n_nd

    end
  end

#corner boundary
  n_count = 1
  for j = 1:(Parameter.N2+1)

  for i = 1:(Parameter.N1+1)
    #festlegen der Boundary: linker Rand= 1; unterer =2; rechter = 3; oberer = 4
    if(i == 1)
      c_boundary[n_count,1] = 1
    end
    if(i==(Parameter.N1+1))
        c_boundary[n_count,3] = 1
    end
    if(j==1)
        c_boundary[n_count,2] = 1
    end
    if(j==(Parameter.N2+1))
        c_boundary[n_count,4] = 1
    end
      n_count += 1
    end
  end
  e_count = 0
  for j = 1:Parameter.N2
    for i = 1:Parameter.N1
      e_count += 1
      corner[e_count,1] = (j-1)* (Parameter.N1+1) + i
      corner[e_count,2] = (j-1)* (Parameter.N1+1) + (i + 1)
      corner[e_count,3] = (j)*(Parameter.N1+1) + (i+1)
      corner[e_count,4] = (j)*(Parameter.N1+1) + i

    end
  end
  tau = zeros(3, 4, n_el)
  data = Data(dim, n_el, n_nd*2, n_corner, x_node, square, n_boundary, x_corner, corner, c_boundary,zeros(2*n_nd),zeros(n_corner), tau)
  return data
end

#Ableitungen der Hutfunktion
function shape(x,y,order)

  if(order==1)

    dN=zeros(2,4)
    dN= [y-1 1-y y -y ; x-1 -x x 1-x]
    NN=zeros(1,4)
    NN=[(1-x)*(1-y) x*(1-y) x*y (1-x)*y ]

    return dN,NN

  elseif(order==2)
    dN=zeros(2,9)
    NN=zeros(2,9)
    x1 = (1-2*x)*(1-x)
    x2 = 4*x*(1-x)
    x3 = x*(2*x-1)
    dx1 = 4*x-3
    dx2 = 4 - 8*x
    dx3 = 4 * x -1

    y1 = (1-2*y)*(1-y)
    y2 = 4*y*(1-y)
    y3 = y*(2*y-1)

    dy1 = 4*y-3
    dy2 = 4 - 8*y
    dy3 = 4 * y -1

    NN=[x1 x3 x3 x1 x2 x3 x2 x1 x2;
        y1 y1 y3 y3 y1 y2 y3 y2 y2]

       dN=[ dx1*y1 dx3*y1 dx3*y3 dx1*y3 dx2*y1 dx3*y2 dx2*y3 dx1*y2 dx2*y2;
              dy1*x1 dy1*x3 dy3*x3 dy3*x1 dy1*x2 dy2*x3 dy3*x2 dy2*x1 dy2*x2]

  return dN,NN

  else
    Println("Not suported order")
  end

end

function gausslegendrerule(n::Integer)
  λ, V = eigen(.5 * SymTridiagonal(ones(n), 1 ./sqrt.(4 .- 1 ./ (1:(n-1)).^2)));
  return λ, V[1,:].^2
end
