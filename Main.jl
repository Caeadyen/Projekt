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
    rhs_vec=zeros(data.n_nd+data.n_corner)
    rh=zeros(data.n_nd+data.n_corner)
    VP=zeros(data.n_nd+data.n_corner)


    for j=1:1
    #Zusammensetzen des Gleichungsystem, dabei wird ignoriert das der Rand 0 ist
        Kvv, Kvp, Kpp, F = assemble_whole(data ,parameter)
        B=[Kvv Kvp'; Kvp Kpp]
        rhs_vec = [F;zeros(data.n_corner)]
        #K=[Kvv Kvp';Kvp Kpp]
        p=zeros(data.n_corner)
        u=zeros(data.n_nd)
        #FF=[F;zeros(data.n_corner)]
        #uu=K\FF
        K=sparse(Kvv+(Kvp'*inv(Matrix(Kpp))*Kvp))
        # for i=1:10
        #     u= K\(F-Kvp'*p)
        #     divu=Kvp*u
        #     println(maximum(divu))
        #     p+=inv(Matrix(Kpp))*divu
        # end
        for i=1:5
            rh= rhs_vec - B*VP
            #rh[data.n_nd+1:end]=rhs_vec[data.n_nd+1:end]-B[data.n_nd+1:end,1:data.n_nd]*VP[1:data.n_nd]
            dVP = B\rh
            VP = VP + dVP
            u=VP[1:data.n_nd]
            p=VP[data.n_nd+1:end]
            div=Kvp*u
            div_max=maximum(abs.(div))
            #println(div_max)
        end


        data.u+=u
        data.p+=p
    end

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
    #Kvv, Kvp, F = assemble_whole(data ,parameter)
    # K = sparse(-Kvp*inv(Matrix(Kvv))*Kvp')
    #     for i=1:10
    #
    #         ru=F
    #         uf=Kvv\ru
    #         pd=K\(rp-Kvp*uf)
    #         #println(pd)
    #         ul=Kvv\(Kvp'*pd)
    #         #println(ul)
    #         ud=uf-ul
    #         u+=ud
    #         p+=pd
    #         ru=F-Kvv*u-Kvp'*p
    #         rp=-Kvp*u
    #
    #         data.u=u
    #         #println(F)
    #         println(norm(rp))
    #         println(norm(ru))
    #         #println(u)
    #
    #     end

    #data.u=VP[1:data.n_nd]
    data.x_node[1:2:end,1]+=data.u[1:2:end]
    data.x_node[2:2:end,1]+=data.u[1:2:end]
    data.x_node[1:2:end,2]+=data.u[2:2:end]
    data.x_node[2:2:end,2]+=data.u[2:2:end]
    #println(norm(rp))
    #println(norm(ru))
    #println(p)
    #println(u)
    #die Punkte des Grids verschieben

    # for i=1:parameter.N2+1
    #     data.x_corner[(1+((i-1)*(parameter.N1+1))):(i*(parameter.N1+1)),2]+=data.u[(1+2*(i-1)*2*(2*parameter.N1+1)):4:(2*(2*parameter.N1+1)+((i-1)*4*(2*parameter.N1+1)))]
    #     data.x_corner[(1+((i-1)*(parameter.N1+1))):(i*(parameter.N1+1)),2]+=data.u[(2+2*(i-1)*2*(2*parameter.N1+1)):4:(2*(2*parameter.N1+1)+((i-1)*4*(2*parameter.N1+1)))]
    # end
    #println(uu)
    return data
end
