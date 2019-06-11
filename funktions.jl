using SparseArrays
using LinearAlgebra
include("tools.jl")
#Building the grid
function builddata(n1,n2)

  h1=1/(2*n1)
  h2=1/(2*n2)
  dim = 2
  n_el = n1*n2
  n_nd = ((2*n1)+1)*((2*n2)+1)*2
  n_corner = (n1+1)*(n2+1)
  n_boundary = zeros(Int,2*n_nd,4)
  x_node = zeros(Float64,2*n_nd,2)
  c_boundary = zeros(Int,n_corner,4)
  x_corner = zeros(Float64,n_corner,2)
  square = zeros(Int,n_el,18)
  corner = zeros(Int,n_el,4)
  n_count = 1

for j = 1:((2*n2)+1)
  yj = (j-1)*h2
  for i = 1:((2*n1)+1)
    xi = (i-1)*h1
    #festlegen der Boundary: linker Rand= 1; unterer =2; rechter = 3; oberer = 4
    x_node[n_count,1] = xi
    x_node[n_count+1,1] = xi
    x_node[n_count,2] = yj
    x_node[n_count+1,2] = yj
    if(i == 1)
      n_boundary[n_count,1] = 1
      n_boundary[n_count+1,1] = 1
    elseif(i==(2*n1)+1)
      n_boundary[n_count,3] = 1
      n_boundary[n_count+1,3] = 1
    elseif(j==1)
      n_boundary[n_count,2] = 1
      n_boundary[n_count+1,2] = 1
    elseif(j==(2*n2)+1)
      n_boundary[n_count,4] = 1
      n_boundary[n_count+1,4] = 1
    end

    n_count += 2
  end
end




e_count = 0

for j = 2:2:(2*n2)
  for i = 1:2:(2*n1-1)
    e_count += 1
    square[e_count,1] = (j-2)* ((4*n1)+2) + (2*i)-1
    square[e_count,2] = (j-2)* ((4*n1)+2) + (2*i)

    square[e_count,3] = (j-2)* ((4*n1)+2) + (2*i)+3
    square[e_count,4] = (j-2)* ((4*n1)+2) + (2*i)+4

    square[e_count,5] = (j)*((4*n1)+2) + (2*i)+3
    square[e_count,6] = (j)*((4*n1)+2) + (2*i)+4

    square[e_count,7] = (j)*((4*n1)+2) + (2*i)-1
    square[e_count,8] = (j)*((4*n1)+2) + (2*i)

    square[e_count,9] = (j-2)* ((4*n1)+2) + (2*i)+1
    square[e_count,10] = (j-2)* ((4*n1)+2) + (2*i)+2

    square[e_count,11] = (j-1)* ((4*n1)+2) + (2*i)+3
    square[e_count,12] = (j-1)* ((4*n1)+2) + (2*i)+4

    square[e_count,13] = (j)* ((4*n1)+2) + (2*i)+1
    square[e_count,14] = (j)* ((4*n1)+2) + (2*i)+2

    square[e_count,15] = (j-1)* ((4*n1)+2) + (2*i)-1
    square[e_count,16] = (j-1)* ((4*n1)+2) + (2*i)

    square[e_count,17] = (j-1)* ((4*n1)+2) + (2*i)+1
    square[e_count,18] = (j-1)* ((4*n1)+2) + (2*i)+2

    end
end

#corner grid
n_count = 1
for j = 1:(n2+1)
  yj = (j-1)*h2*2
  for i = 1:(n1+1)
    xi = (i-1)*h1*2
    #festlegen der Boundary: linker Rand= 1; unterer =2; rechter = 3; oberer = 4
    x_corner[n_count,1] = xi
    x_corner[n_count,2] = yj

    if(i == 1)
      c_boundary[n_count,1] = 1

    elseif(i==(n1+1))
      c_boundary[n_count,3] = 1

    elseif(j==1)
      c_boundary[n_count,2] = 1

    elseif(j==(n2+1))
      c_boundary[n_count,4] = 1

    end

    n_count += 1
  end
end
e_count = 0
for j = 1:n2
  for i = 1:n1
    e_count += 1
    corner[e_count,1] = (j-1)* (n1+1) + i
    corner[e_count,2] = (j-1)* (n1+1) + (i + 1)
    corner[e_count,3] = (j)*(n1+1) + (i+1)
    corner[e_count,4] = (j)*(n1+1) + i

    end
end
#u = zeros(2*data.n_nd)
#  p = zeros(data.n_corner)
data = Data(dim, n_el, n_nd, n_corner, x_node, square, n_boundary, x_corner, corner, c_boundary,zeros(2*n_nd),zeros(n_corner))
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
  #NN=zeros(2,18)
    NN=[x1*y1 0 x3*y1 0 x3*y3 0 x1*y3 0 x2*y1 0 x3*y2 0 x2*y3 0 x1*y2 0 x2*y2 0;
       0 x1*y1 0 x3*y1 0 x3*y3 0 x1*y3 0 x2*y1 0 x3*y2 0 x2*y3 0 x1*y2 0 x2*y2]

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
