include("funktions.jl")
include("tools.jl")
include("solver.jl")
function solvestoke(parameter)

    data = builddata(parameter.N1,parameter.N2)
    #in data machen
    #TODO versuch, muss noch ne schleife werden
    #p = zeros(data.n_corner)
    uu=[data.u;data.p]

    #Zusammensetzen des Gleichungsystem, dabei wird ignoriert das der Rand 0 ist
    Kvv, Kvp, Kpp, F = assemble_whole(data ,parameter)

    K=[Kvv Kvp';Kvp -Kpp]

    FF=[F;-Kpp*data.p]

    #Finden der inneren Knoten, des Displacment teils
    inner_node=Int64[]
    for j = 1:(data.n_nd)

        if(parameter.boundary*data.n_boundary[j,:]==[0])

            push!(inner_node,j)

        end

    end
    #Finden der inneren Knoten, des Druck teils
    for j = 1:data.n_corner

        if(parameter.boundary*data.c_boundary[j,:]==[0])

            push!(inner_node,j+(data.n_nd))

        end

    end

    #Lösen des Gleichungsystems auf den inneren Knoten und ploten des Ergebnisses

    uu[inner_node] = K[inner_node,inner_node]\FF[inner_node]

    
    #TODO verschiebung der punkte
    u=uu[1:data.n_nd]
    p=uu[data.n_nd+1:end]


end
