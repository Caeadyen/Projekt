using PyPlot
include("funktions.jl")
include("tools.jl")
include("solver.jl")
function solvestoke(parameter)

    data = builddata(parameter.N1,parameter.N2)
    #in data machen
    #TODO versuch, muss noch ne schleife werden
    p = zeros(data.n_corner)
    uu=zeros(length(data.n_nd))
    for i=1:1
    #Zusammensetzen des Gleichungsystem, dabei wird ignoriert das der Rand 0 ist
        Kvv, Kvp, Kpp, F = assemble_whole(data ,parameter)
        #K=[Kvv Kvp';Kvp -Kpp]

        #FF=[F;-Kpp*data.p]

    #Finden der inneren Knoten, des Displacment teils
        inner_node=Int64[]
        for j = 1:(data.n_nd)
            if(data.n_boundary[j,2]==0&&data.n_boundary[j,4]==0)

                push!(inner_node,j)

            end

        end
    #Finden der inneren Knoten, des Druck teils
        for j = 1:data.n_corner

            if(data.c_boundary[j,2]==0&&data.c_boundary[j,4]==0)

                push!(inner_node,j+(data.n_nd))

            end

        end

    #Lösen des Gleichungsystems auf den inneren Knoten und ploten des Ergebnisses

    #uu[inner_node] = K[inner_node,inner_node]\FF[inner_node]
        K=sparse(Kvv+(Kvp'*inv(Matrix(Kpp))*Kvp))
        uu=K\(F-Kvp'*p)
        #println(uu)
        divv=Kvp*uu
        p=Kpp\(Kpp*p+divv)
        data.u+=uu

    end
    #die Punkte des Grids verschieben
    data.x_node[1:2:end,1]+=data.u[1:2:end]
    data.x_node[2:2:end,1]+=data.u[1:2:end]
    data.x_node[1:2:end,2]+=data.u[2:2:end]
    data.x_node[2:2:end,2]+=data.u[2:2:end]
    for i=1:parameter.N2+1
        data.x_corner[(1+((i-1)*(parameter.N1+1))):(i*(parameter.N1+1)),2]+=data.u[(1+2*(i-1)*2*(2*parameter.N1+1)):4:(2*(2*parameter.N1+1)+((i-1)*4*(2*parameter.N1+1)))]
        data.x_corner[(1+((i-1)*(parameter.N1+1))):(i*(parameter.N1+1)),2]+=data.u[(2+2*(i-1)*2*(2*parameter.N1+1)):4:(2*(2*parameter.N1+1)+((i-1)*4*(2*parameter.N1+1)))]
    end
    #println(uu)
    return data
end
