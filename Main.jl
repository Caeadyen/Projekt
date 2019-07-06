using PyPlot
include("funktions.jl")
include("tools.jl")
include("solver.jl")
function solvestoke(parameter)

    data = builddata(parameter)

    for j=1:parameter.dt:parameter.maxtime
    #Zusammensetzen des Gleichungsystem, dabei wird ignoriert das der Rand 0 ist
        Kvv, Kvp, F = assemble_whole(data ,parameter)

    #Finden der inneren Knoten, des Displacment teils
      inner_node=Int64[]
      boundary_nodel=Int64[]
      boundary_noder=Int64[]

       for j = 1:(data.n_nd)
           if(data.n_boundary[j,1]==0&&data.n_boundary[j,3]==0)

               push!(inner_node,j)
           else
               if(data.n_boundary[j,1]==1)
                   push!(boundary_nodel,j)
               else
                   push!(boundary_noder,j)
               end
           end

       end
       yeartime=3600*24*365.25
       FF=[F;zeros(data.n_corner)]
       K=[Kvv Kvp';Kvp spzeros(data.n_corner,data.n_corner)]
       boundary_valuel=parameter.lboundary*ones(size(boundary_nodel,1))/yeartime
       boundary_valuer=parameter.rboundary*ones(size(boundary_noder,1))/yeartime

       #Applying Boundray, just pushing from both sides
       FF=FF-K[:,boundary_nodel]*boundary_valuel;
       FF=FF-K[:,boundary_noder]*boundary_valuer;
       uu=zeros(data.n_nd+data.n_corner)

       #solving the system
       uu[inner_node] = K[inner_node,inner_node]\FF[inner_node]
       uu[boundary_nodel] = boundary_valuel
       uu[boundary_noder] = boundary_valuer

        dt=3600*24*365.25 #100jahre
        data.u[inner_node]+=uu[1:size(inner_node,1)]*parameter.dt*yeartime
        data.x_node[1:end,1]+=uu[1:Int(data.n_nd/2)]*parameter.dt*yeartime
        data.x_node[1:end,2]+=uu[Int(data.n_nd/2)+1:data.n_nd]*parameter.dt*yeartime

    end

    return data
end
