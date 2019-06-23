using PyPlot
include("funktions.jl")
include("tools.jl")
include("solver.jl")
function solvestoke(parameter)

    data = builddata(parameter.N1,parameter.N2)
    #in data machen
    #TODO versuch, muss noch ne schleife werden
    rp = zeros(data.n_corner)
    rd = zeros(data.n_corner)
    ru=zeros(data.n_nd)
    uf=zeros(data.n_nd)
    ul=zeros(data.n_nd)
    u=zeros(data.n_nd)
    p=zeros(data.n_corner)

    #Zusammensetzen des Gleichungsystem, dabei wird ignoriert das der Rand 0 ist
        #Kvv, Kvp, F = assemble_whole(data ,parameter)

        # K=[Kvv Kvp';Kvp zeros(data.n_corner,data.n_corner)]
        #
        # FF=[F;zeros(data.n_corner)]
        # uu=K\FF
        # u=uu[1:data.n_nd]
    #Finden der inneren Knoten, des Displacment teils
    #    inner_node=Int64[]
    #    for j = 1:(data.n_nd)
    #        if(data.n_boundary[j,2]==0&&data.n_boundary[j,4]==0)

    #            push!(inner_node,j)

    #        end

    #    end
    #Finden der inneren Knoten, des Druck teils
    #    for j = 1:data.n_corner

    #        if(data.c_boundary[j,2]==0&&data.c_boundary[j,4]==0)

    #            push!(inner_node,j+(data.n_nd))

    #        end

    #    end

    #Lösen des Gleichungsystems auf den inneren Knoten und ploten des Ergebnisses
    #println(F)
    #uu[inner_node] = K[inner_node,inner_node]\FF[inner_node]

    for i=1:10
        Kvv, Kvp, F = assemble_whole(data ,parameter)
        K = sparse(-Kvp*inv(Matrix(Kvv))*Kvp')
        ru=F
        uf=Kvv\ru
        pd=K\(rp-Kvp*uf)
        #println(pd)
        ul=Kvv\(Kvp'*pd)
        #println(ul)
        ud=uf-ul
        u+=ud
        p+=pd
        ru=F-Kvv*u-Kvp'*p
        rp=-Kvp*u

        data.u=u
        #println(F)
        println(norm(rp))
        println(norm(ru))
        #println(u)

    end
    data.x_node[1:2:end,1]+=u[1:2:end]
    data.x_node[2:2:end,1]+=u[1:2:end]
    data.x_node[1:2:end,2]+=u[2:2:end]
    data.x_node[2:2:end,2]+=u[2:2:end]
    #println(norm(rp))
    #println(norm(ru))
    #println(p)
    #println(u)
    #die Punkte des Grids verschieben

    for i=1:parameter.N2+1
        data.x_corner[(1+((i-1)*(parameter.N1+1))):(i*(parameter.N1+1)),2]+=data.u[(1+2*(i-1)*2*(2*parameter.N1+1)):4:(2*(2*parameter.N1+1)+((i-1)*4*(2*parameter.N1+1)))]
        data.x_corner[(1+((i-1)*(parameter.N1+1))):(i*(parameter.N1+1)),2]+=data.u[(2+2*(i-1)*2*(2*parameter.N1+1)):4:(2*(2*parameter.N1+1)+((i-1)*4*(2*parameter.N1+1)))]
    end
    #println(uu)
    return data
end
