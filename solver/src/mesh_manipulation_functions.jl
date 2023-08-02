using SparseArrays
#Aggiornata
function From_3D_to_1D(i, j, k, M, N) 
	pos = ((k-1) * M * N) + ((j-1) * M) + i;
    return convert(Int64,pos)
end

#Aggiornata
function bin_search(num, A)
    index = 0
    n = length(A)
    left = 1
    right = n
    while left <= right
        mid = ceil((left + right) / 2)
        if A[mid] == num
            index = mid
            break
        else
            if A[mid] > num
                right = mid - 1
            else
                left = mid + 1
            end
        end
    end
    return index
end

#Non sostituita perché identica alle corrispondenti
function create_volumes_mapping_and_centers(matrice,Nx,Ny,Nz,num_centri,sx,sy,sz,min_v)


    println("----",(Nx * Ny * Nz))
 
    mapping = zeros(Int64, Nx * Ny * Nz)
    centri_vox = zeros(Float64, num_centri, 3)
    id_mat = zeros(Int64, num_centri)

    num_grids = size(matrice)[1]
    num_ele=1
    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if (matrice[k][cont][cont2][cont3]==true)
                        mapping[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]=num_ele
                        centri_vox[num_ele, 1] = min_v[1] + (sx * (cont - 1))  + (sx / 2.0)
                        centri_vox[num_ele, 2] = min_v[2] + (sy * (cont2 - 1)) + (sy / 2.0)
                        centri_vox[num_ele, 3] = min_v[3] + (sz * (cont3 - 1)) + (sz / 2.0)
                        id_mat[num_ele] = k
                        num_ele = num_ele + 1
                        break
                    end
                end
            end
        end
    end


    return num_ele,mapping,centri_vox,id_mat
end

#Aggiornata
function create_nodes_ref(grids, Nx, Ny, Nz, num_full_vox, external_grids, mapping_vols, dominant_list)
    num_grids = size(external_grids, 1)
    nodes = zeros(num_full_vox)
    lista_nod = zeros(num_full_vox)
    cont_d = 0
    for ka in range(1,length(dominant_list))
        k = dominant_list[ka]
        for cont = 1:Nx
            for cont2 = 1:Ny
                for cont3 = 1:Nz
                    if grids[k][cont][cont2][cont3]
                        c1 = 1 + 2 * (cont - 1) + 1
                        c2 = 1 + 2 * (cont2 - 1) + 1
                        c3 = 1 + 2 * (cont3 - 1) + 1
                        f1 = external_grids[k,1][cont,cont2,cont3]
                        f2 = external_grids[k,2][cont,cont2,cont3]
                        f3 = external_grids[k,3][cont,cont2,cont3]
                        f4 = external_grids[k,4][cont,cont2,cont3]
                        f5 = external_grids[k,5][cont,cont2,cont3]
                        f6 = external_grids[k,6][cont,cont2,cont3]
                        f1_c = false
                        f2_c = false
                        f3_c = false
                        f4_c = false
                        f5_c = false
                        f6_c = false
                        if cont2 - 1 > 1
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,2][cont,cont2 - 1,cont3]
                                        f1_c = true
                                    end
                                end
                            end
                        end
                        if cont2 + 1 <= Ny
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,1][cont,cont2 + 1,cont3]
                                        f2_c = true
                                    end
                                end
                            end
                        end
                        if cont - 1 > 1
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,4][cont - 1,cont2,cont3]
                                        f3_c = true
                                    end
                                end
                            end
                        end
                        if cont + 1 <= Nx
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,3][cont + 1,cont2,cont3]
                                        f4_c = true
                                    end
                                end
                            end
                        end
                        if cont3 - 1 > 1
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,6][cont,cont2,cont3 - 1]
                                        f5_c = true
                                    end
                                end
                            end
                        end
                        if cont3 + 1 <= Nz
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,5][cont,cont2,cont3 + 1]
                                        f6_c = true
                                    end
                                end
                            end
                        end
                        is_f1 = f1
                        is_f2 = f2
                        if f1 && f2
                            is_f1 = false
                            is_f2 = false
                            if f1_c && !f2_c
                                is_f1 = true
                            elseif !f1_c && f2_c
                                is_f2 = true
                            end
                        end
                        is_f3 = f3
                        is_f4 = f4
                        if f3 && f4
                            is_f3 = false
                            is_f4 = false
                            if f3_c && !f4_c
                                is_f3 = true
                            elseif !f3_c && f4_c
                                is_f4 = true
                            end
                        end
                        is_f5 = f5
                        is_f6 = f6
                        if f5 && f6
                            is_f5 = false
                            is_f6 = false
                            if f5_c && !f6_c
                                is_f5 = true
                            elseif !f5_c && f6_c
                                is_f6 = true
                            end
                        end
                        if any([is_f1 is_f2 is_f3 is_f4 is_f5 is_f6])
                            if is_f1 && !is_f2 && !is_f3 && !is_f4 && !is_f5 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2 - 1, c3, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f2 && !is_f1 && !is_f3 && !is_f4 && !is_f5 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2 + 1, c3, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f3 && !is_f1 && !is_f2 && !is_f4 && !is_f5 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 - 1, c2, c3, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f4 && !is_f1 && !is_f2 && !is_f3 && !is_f5 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 + 1, c2, c3, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f5 && !is_f1 && !is_f2 && !is_f3 && !is_f4 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2, c3 - 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f6 && !is_f1 && !is_f2 && !is_f3 && !is_f4 && !is_f5
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2, c3 + 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f1 && is_f3 && !is_f2 && !is_f4 && !is_f5 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 - 1, c2 - 1, c3, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f1 && is_f4 && !is_f2 && !is_f3 && !is_f5 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 + 1, c2 - 1, c3, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f1 && is_f5 && !is_f2 && !is_f3 && !is_f4 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2 - 1, c3 - 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f1 && is_f6 && !is_f2 && !is_f3 && !is_f4 && !is_f5
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2 - 1, c3 + 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f1 && is_f3 && is_f5 && !is_f2 && !is_f4 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 - 1, c2 - 1, c3 - 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f1 && is_f3 && is_f6 && !is_f2 && !is_f4 && !is_f5
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 - 1, c2 - 1, c3 + 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f1 && is_f4 && is_f5 && !is_f2 && !is_f3 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 + 1, c2 - 1, c3 - 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f1 && is_f4 && is_f6 && !is_f2 && !is_f3 && !is_f5
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 + 1, c2 - 1, c3 + 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f2 && is_f3 && !is_f1 && !is_f4 && !is_f5 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 - 1, c2 + 1, c3, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f2 && is_f4 && !is_f1 && !is_f3 && !is_f5 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 + 1, c2 + 1, c3, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f2 && is_f5 && !is_f1 && !is_f3 && !is_f4 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2 + 1, c3 - 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f2 && is_f6 && !is_f1 && !is_f3 && !is_f4 && !is_f5
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2 + 1, c3 + 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f2 && is_f3 && is_f5 && !is_f1 && !is_f4 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 - 1, c2 + 1, c3 - 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f2 && is_f3 && is_f6 && !is_f1 && !is_f4 && !is_f5
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 - 1, c2 + 1, c3 + 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f2 && is_f4 && is_f5 && !is_f1 && !is_f3 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 + 1, c2 + 1, c3 - 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f2 && is_f4 && is_f6 && !is_f1 && !is_f3 && !is_f5
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 + 1, c2 + 1, c3 + 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f3 && is_f5 && !is_f1 && !is_f2 && !is_f4 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 - 1, c2, c3 - 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f3 && is_f6 && !is_f1 && !is_f2 && !is_f4 && !is_f5
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 - 1, c2, c3 + 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f4 && is_f5 && !is_f1 && !is_f2 && !is_f3 && !is_f6
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 + 1, c2, c3 - 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            elseif is_f4 && is_f6 && !is_f1 && !is_f2 && !is_f3 && !is_f5
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1 + 1, c2, c3 + 1, 3 * Nx, 3 * Ny)
                                cont_d += 1
                                lista_nod[cont_d] = nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])]
                            end
                        else
                            nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2, c3, 3 * Nx, 3 * Ny)
                        end
                    end
                end
            end
        end
    end
    nodes_red_dom = sort(unique(lista_nod[1:cont_d]))
    non_domin_list = setdiff(1:num_grids, dominant_list)
    nodes_reused = zeros(num_full_vox)
    cont_reu = 0
    for ka in range(1,length(non_domin_list))
        k = non_domin_list[ka]
        for cont = 1:Nx
            for cont2 = 1:Ny
                for cont3 = 1:Nz
                    if grids[k][cont][cont2][cont3]
                        c1 = 1 + 2 * (cont - 1) + 1
                        c2 = 1 + 2 * (cont2 - 1) + 1
                        c3 = 1 + 2 * (cont3 - 1) + 1
                        f1 = external_grids[k,1][cont,cont2,cont3]
                        f2 = external_grids[k,2][cont,cont2,cont3]
                        f3 = external_grids[k,3][cont,cont2,cont3]
                        f4 = external_grids[k,4][cont,cont2,cont3]
                        f5 = external_grids[k,5][cont,cont2,cont3]
                        f6 = external_grids[k,6][cont,cont2,cont3]
                        f1_c = false
                        f2_c = false
                        f3_c = false
                        f4_c = false
                        f5_c = false
                        f6_c = false
                        if cont2 - 1 > 1
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,2][cont,cont2 - 1,cont3]
                                        f1_c = true
                                    end
                                end
                            end
                        end
                        if cont2 + 1 <= Ny
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,1][cont,cont2 + 1,cont3]
                                        f2_c = true
                                    end
                                end
                            end
                        end
                        if cont - 1 > 1
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,4][cont - 1,cont2,cont3]
                                        f3_c = true
                                    end
                                end
                            end
                        end
                        if cont + 1 <= Nx
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,3][cont + 1,cont2,cont3]
                                        f4_c = true
                                    end
                                end
                            end
                        end
                        if cont3 - 1 > 1
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,6][cont,cont2,cont3 - 1]
                                        f5_c = true
                                    end
                                end
                            end
                        end
                        if cont3 + 1 <= Nz
                            for k2 = 1:num_grids
                                if k != k2
                                    if external_grids[k2,5][cont,cont2,cont3 + 1]
                                        f6_c = true
                                    end
                                end
                            end
                        end
                        is_f1 = f1
                        is_f2 = f2
                        if f1 && f2
                            is_f1 = false
                            is_f2 = false
                            if f1_c && !f2_c
                                is_f1 = true
                            elseif !f1_c && f2_c
                                is_f2 = true
                            end
                        end
                        is_f3 = f3
                        is_f4 = f4
                        if f3 && f4
                            is_f3 = false
                            is_f4 = false
                            if f3_c && !f4_c
                                is_f3 = true
                            elseif !f3_c && f4_c
                                is_f4 = true
                            end
                        end
                        is_f5 = f5
                        is_f6 = f6
                        if f5 && f6
                            is_f5 = false
                            is_f6 = false
                            if f5_c && !f6_c
                                is_f5 = true
                            elseif !f5_c && f6_c
                                is_f6 = true
                            end
                        end
                        if any([is_f1, is_f2, is_f3, is_f4, is_f5, is_f6])
                            nodes_to_see = build_nodes(c1, c2, c3, 3*Nx, 3*Ny)
                            nodo_shared, val_nodo = bin_search_mod(nodes_to_see, nodes_red_dom)
                            if abs(nodo_shared) > 1e-8
                                nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = val_nodo
                                cont_reu += 1
                                nodes_reused[cont_reu] = val_nodo
                            else
                                if is_f1 && !is_f2 && !is_f3 && !is_f4 && !is_f5 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2-1, c3, 3*Nx, 3*Ny)
                                elseif is_f2 && !is_f1 && !is_f3 && !is_f4 && !is_f5 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2+1, c3, 3*Nx, 3*Ny)
                                elseif is_f3 && !is_f1 && !is_f2 && !is_f4 && !is_f5 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1-1, c2, c3, 3*Nx, 3*Ny)
                                elseif is_f4 && !is_f1 && !is_f2 && !is_f3 && !is_f5 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1+1, c2, c3, 3*Nx, 3*Ny)
                                elseif is_f5 && !is_f1 && !is_f2 && !is_f3 && !is_f4 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2, c3-1, 3*Nx, 3*Ny)
                                elseif is_f6 && !is_f1 && !is_f2 && !is_f3 && !is_f4 && !is_f5
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2, c3+1, 3*Nx, 3*Ny)
                                elseif is_f1 && is_f3 && !is_f2 && !is_f4 && !is_f5 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1-1, c2-1, c3, 3*Nx, 3*Ny)
                                elseif is_f1 && is_f4 && !is_f2 && !is_f3 && !is_f5 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1+1, c2-1, c3, 3*Nx, 3*Ny)
                                elseif is_f1 && is_f5 && !is_f2 && !is_f3 && !is_f4 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2-1, c3-1, 3*Nx, 3*Ny)
                                elseif is_f1 && is_f6 && !is_f2 && !is_f3 && !is_f4 && !is_f5
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2-1, c3+1, 3*Nx, 3*Ny)
                                elseif is_f1 && is_f3 && is_f5 && !is_f2 && !is_f4 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1-1, c2-1, c3-1, 3*Nx, 3*Ny)
                                elseif is_f1 && is_f3 && is_f6 && !is_f2 && !is_f4 && !is_f5
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1-1, c2-1, c3+1, 3*Nx, 3*Ny)
                                elseif is_f1 && is_f4 && is_f5 && !is_f2 && !is_f3 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1+1, c2-1, c3-1, 3*Nx, 3*Ny)
                                elseif is_f1 && is_f4 && is_f6 && !is_f2 && !is_f3 && !is_f5
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1+1, c2-1, c3+1, 3*Nx, 3*Ny)
                                elseif is_f2 && is_f3 && !is_f1 && !is_f4 && !is_f5 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1-1, c2+1, c3, 3*Nx, 3*Ny)
                                elseif is_f2 && is_f4 && !is_f1 && !is_f3 && !is_f5 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1+1, c2+1, c3, 3*Nx, 3*Ny)
                                elseif is_f2 && is_f5 && !is_f1 && !is_f3 && !is_f4 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2+1, c3-1, 3*Nx, 3*Ny)
                                elseif is_f2 && is_f6 && !is_f1 && !is_f3 && !is_f4 && !is_f5
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2+1, c3+1, 3*Nx, 3*Ny)
                                elseif is_f2 && is_f3 && is_f5 && !is_f1 && !is_f4 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1-1, c2+1, c3-1, 3*Nx, 3*Ny)
                                elseif is_f2 && is_f3 && is_f6 && !is_f1 && !is_f4 && !is_f5
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1-1, c2+1, c3+1, 3*Nx, 3*Ny)
                                elseif is_f2 && is_f4 && is_f5 && !is_f1 && !is_f3 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1+1, c2+1, c3-1, 3*Nx, 3*Ny)
                                elseif is_f2 && is_f4 && is_f6 && !is_f1 && !is_f3 && !is_f5
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1+1, c2+1, c3+1, 3*Nx, 3*Ny)
                                elseif is_f3 && is_f5 && !is_f1 && !is_f2 && !is_f4 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1-1, c2, c3-1, 3*Nx, 3*Ny)
                                elseif is_f3 && is_f6 && !is_f1 && !is_f2 && !is_f4 && !is_f5
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1-1, c2, c3+1, 3*Nx, 3*Ny)
                                elseif is_f4 && is_f5 && !is_f1 && !is_f2 && !is_f3 && !is_f6
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1+1, c2, c3-1, 3*Nx, 3*Ny)
                                elseif is_f4 && is_f6 && !is_f1 && !is_f2 && !is_f3 && !is_f5
                                    nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1+1, c2, c3+1, 3*Nx, 3*Ny)
                                end
                            end
                        else
                            nodes[convert(Int64,mapping_vols[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])] = From_3D_to_1D(c1, c2, c3, 3*Nx, 3*Ny)
                        end
                    end
                end
            end
        end
    end

    nodes_red = sort(unique(nodes))
    nodes_reused_clean = sort(unique(nodes_reused[1:cont_reu]))

    return nodes,nodes_red,nodes_reused_clean
end

#Aggiornata
function distfcm(center, data)
    out = zeros(size(center, 1), size(data, 1))
    if size(center[1][1], 2) > 1
        for k in range(1,size(center[1][1], 1))
            out[k, :] = sqrt.(sum((data .- ones(size(data, 1), 1) .* center[k][k]).^2, dims=2))
        end
    else
        for k in range(1,size(center, 1))
            out[k, :] = transpose(abs.(center[k] .- data))
        end
    end
    return out
end

#Aggiornata
function nodes_find_rev(Nodes_inp_coord, nodi_centri, node_to_skip)
    indici = sortperm(vec(distfcm(Nodes_inp_coord, nodi_centri)))
    if indici[1] != node_to_skip
        nodes = indici[1]
    else
        nodes = indici[2]
    end
    return nodes
end

#Aggiornata
function find_nodes_port(nodi_centri, port_start, port_end, nodi, nodi_red)
    N = size(port_start)[1]
    port_voxels = zeros(Int64, N, 2)
    port_nodes = zeros(Int64, N, 2)
    for cont in range(1, stop=N)
        port_voxels[cont, 1] = nodes_find_rev(port_start[cont,:], nodi_centri, -1)
        port_voxels[cont, 2] = nodes_find_rev(port_end[cont,:],nodi_centri, port_voxels[cont,1])
        port_nodes[cont, 1] = bin_search(nodi[convert(Int64,port_voxels[cont, 1])], nodi_red)
        port_nodes[cont, 2] = bin_search(nodi[convert(Int64,port_voxels[cont, 2])], nodi_red)
    end
    return port_voxels, port_nodes
end

#Non sostituita
function create_external_grids(matrice,Nx,Ny,Nz)
    num_grids = size(matrice)[1]
    
    OUTPUTgrids = zeros(Int8 ,num_grids ,6,Nx ,Ny ,Nz)


    for cont2 in range(1, stop=Ny)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if (matrice[k][1][cont2][cont3]==1)
                    OUTPUTgrids[k,3,1,cont2,cont3]=1
                end
            end
        end
    end

    for cont2 in range(1, stop=Ny)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if (matrice[k][Nx][cont2][cont3]==1)
                    OUTPUTgrids[k,4,Nx,cont2,cont3]=1
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if (matrice[k][cont][1][cont3]==1)
                    OUTPUTgrids[k,1,cont,1,cont3] = 1
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if (matrice[k][cont][Ny][cont3] == 1)
                    OUTPUTgrids[k,2,cont,Ny,cont3] = 1
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for k in range(1, stop=num_grids)
                if (matrice[k][cont][cont2][1]==1)
                    OUTPUTgrids[k,5,cont,cont2,1] = 1
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for k in range(1, stop=num_grids)
                if(matrice[k][cont][cont2][Nz]==1)
                    OUTPUTgrids[k,6,cont,cont2,Nz] = 1
                end
            end
        end
    end

    for cont in range(2,stop=Nx-1)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if (matrice[k][cont][cont2][cont3] == 1)
                        if (matrice[k][cont - 1][cont2][cont3] == 0)
                            OUTPUTgrids[k,3,cont,cont2,cont3] = 1
                        end
                        if (matrice[k][cont + 1][cont2][cont3] == 0)
                            OUTPUTgrids[k,4,cont,cont2,cont3] = 1
                        end
                    end
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(2, stop=Ny-1)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if (matrice[k][cont][cont2][cont3] == 1)
                        if (matrice[k][cont][cont2 - 1][cont3] == 0)
                            OUTPUTgrids[k,1,cont,cont2,cont3] = 1
                        end
                        if (matrice[k][cont][cont2 + 1][cont3] == 0)
                            OUTPUTgrids[k,2,cont,cont2,cont3] = 1
                        end
                    end
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(2, stop=Nz-1)
                for k in range(1, stop=num_grids)
                    if (matrice[k][cont][cont2][cont3] == 1)
                        if (matrice[k][cont][cont2][cont3 - 1] == 0)
                            OUTPUTgrids[k,5,cont,cont2,cont3] = 1
                        end
                        if (matrice[k][cont][cont2][cont3 + 1] == 0)
                            OUTPUTgrids[k,6,cont,cont2,cont3] = 1
                        end
                    end
                end
            end
        end
    end

    return OUTPUTgrids
end

#Aggiornata
function create_mapping_Ax(grids, mapping_Vox, nodes, nodes_red)
    num_grids = length(grids)
    Nx = size(grids[1],1)
    Ny = size(grids[1][1],1)
    Nz = size(grids[1][1][1],1)
    N_max = (Nx-1)*Ny*Nz
    mapping = zeros(Int64, N_max)
    num_ele = 0
    for cont2=1:Ny
        for cont3=1:Nz
            for cont=1:Nx-1
                for k in 1:num_grids
                    if grids[k][cont][cont2][cont3] && grids[k][cont+1][cont2][cont3]
                        nn1 = bin_search(nodes[convert(Int64,mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])], nodes_red)
                        nn2 = bin_search(nodes[convert(Int64,mapping_Vox[From_3D_to_1D(cont+1, cont2, cont3, Nx, Ny)])], nodes_red)
                        if abs(nn1-nn2) > 1e-8
                            kkey = From_3D_to_1D(cont, cont2, cont3, Nx-1, Ny)
                            if mapping[kkey] == 0
                                num_ele += 1
                                mapping[kkey] = num_ele
                            end
                            break
                        end
                    end
                end
            end
        end
    end
    return mapping, num_ele
end

#Aggiornata
function create_mapping_Ay(grids, mapping_Vox, nodes, nodes_red)
    num_grids = length(grids)
    Nx = size(grids[1],1)
    Ny = size(grids[1][1],1)
    Nz = size(grids[1][1][1],1)
    N_max = Nx * (Ny - 1) * Nz
    mapping = zeros(Int64, N_max)
    num_ele = 0
    for cont3 in 1:Nz
        for cont in 1:Nx
            for cont2 in 1:Ny-1
                for k in 1:num_grids
                    if grids[k][cont][cont2][cont3] && grids[k][cont][cont2+1][cont3]
                        nn1 = bin_search(nodes[convert(Int64,mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])], nodes_red)
                        nn2 = bin_search(nodes[convert(Int64,mapping_Vox[From_3D_to_1D(cont, cont2+1, cont3, Nx, Ny)])], nodes_red)
                        if abs(nn1 - nn2) > 1e-8
                            kkey = From_3D_to_1D(cont, cont2, cont3, Nx, Ny-1)
                            if mapping[kkey] == 0
                                num_ele += 1
                                mapping[kkey] = num_ele
                            end
                            break
                        end
                    end
                end
            end
        end
    end
    return mapping, num_ele
end

#Aggiornata
function create_mapping_Az(grids, mapping_Vox, nodes, nodes_red)
    num_grids = length(grids)
    Nx = size(grids[1],1)
    Ny = size(grids[1][1],1)
    Nz = size(grids[1][1][1],1)
    N_max = Nx * Ny * (Nz - 1)
    mapping = zeros(Int64, N_max)
    num_ele = 0
    for cont = 1:Nx
        for cont2 = 1:Ny
            for cont3 = 1:Nz-1
                for k = 1:num_grids
                    if grids[k][cont][cont2][cont3] && grids[k][cont][cont2][cont3+1]
                        nn1 = bin_search(nodes[convert(Int64,mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)])], nodes_red)
                        nn2 = bin_search(nodes[convert(Int64,mapping_Vox[From_3D_to_1D(cont, cont2, cont3+1, Nx, Ny)])], nodes_red)
                        if abs(nn1-nn2) > 1e-8
                            kkey = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                            if mapping[kkey] == 0
                                num_ele += 1
                                mapping[kkey] = num_ele
                            end
                            break
                        end
                    end
                end
            end
        end
    end
    return mapping, num_ele
end

#Non sostituita
function create_A_mats_volInd(matrice,Nx,Ny,Nz,mapping_Vox,mapAx, NAx, mapAy, NAy, mapAz, NAz, sx,sy,sz,min_v,nodi,nodi_red)
 
    num_grids = size(matrice)[1]
    lix_mat = zeros(Int8, NAx, 2)
    lix_border = zeros(Int8, NAx, 2)
    ind_row = zeros(Int64, 2*NAx+2*NAy+2*NAz)
    ind_col = zeros(Int64, 2*NAx+2*NAy+2*NAz)
    vals_A = zeros(Float64, 2*NAx+2*NAy+2*NAz)
    bars_Lp_x = zeros(Float64, NAx, 6)

    num_ele = 0

    for cont2 in range(1, stop=Ny)
        for cont3 in range(1, stop=Nz)
            for cont in range(1, stop=Nx - 1)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont+1][cont2][cont3] == 1))

                        pos = mapAx[From_3D_to_1D(cont, cont2, cont3, Nx - 1, Ny)] + 1
                        
                        ind_row[num_ele+1] = pos
                        ind_col[num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]], nodi_red)

                        bars_Lp_x[pos,1] = min_v[1] + sx * (cont - 1) + sx/2
                        bars_Lp_x[pos,2] = min_v[2] + sy * (cont2 - 1)
                        bars_Lp_x[pos,3] = min_v[3] + sz * (cont3 - 1)

                        vals_A[num_ele+1] = -1.0
                        num_ele = num_ele + 1

                        ind_row[num_ele+1] = pos
                        ind_col[num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont+1, cont2, cont3, Nx, Ny)]], nodi_red)
                        vals_A[num_ele+1] = 1.0
                        num_ele = num_ele + 1

                        bars_Lp_x[pos,4] = bars_Lp_x[pos,1] + sx
                        bars_Lp_x[pos,5] = bars_Lp_x[pos,2] + sy
                        bars_Lp_x[pos,6] = bars_Lp_x[pos,3] + sz

                        lix_mat[pos, 1] = k+1
                        lix_mat[pos, 2] = k+1

                        if cont > 1
                           if (matrice[k][cont-1][cont2][cont3]==0)
                               lix_border[pos, 1] = k + 1
                               bars_Lp_x[pos,1] = bars_Lp_x[pos,1] - sx/2.0
                           end
                        else
                            lix_border[pos, 1] = k + 1
                            bars_Lp_x[pos,1] = bars_Lp_x[pos,1] - sx / 2.0
                        end

                        if cont + 1 == Nx
                            lix_border[pos, 2] = k + 1
                            bars_Lp_x[pos,4] = bars_Lp_x[pos,4] + sx/2.0
                        else
                            if cont + 2 <= Nx
                                if (matrice[k][cont+2][cont2][cont3]==0)
                                    lix_border[pos, 2] = k + 1
                                    bars_Lp_x[pos,4] = bars_Lp_x[pos,4] + sx/2.0
                                end
                            end
                        end
                        break
                    end
                end
            end
        end
    end


    liy_mat = zeros(Int8, NAy, 2)
    liy_border = zeros(Int8, NAy, 2)
    bars_Lp_y = zeros(Float64, NAy, 6)

    starter=num_ele
    num_ele = 0

    for cont3 in range(1, stop=Nz)
        for cont in range(1, stop=Nx)
            for cont2 in range(1, stop=Ny - 1)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont][cont2+1][cont3] == 1))

                        pos = mapAy[From_3D_to_1D(cont, cont2, cont3, Nx, Ny - 1)] + 1

                        bars_Lp_y[pos,1] = min_v[1] + sx * (cont - 1)
                        bars_Lp_y[pos,2] = min_v[2] + sy * (cont2 - 1) + sy/2
                        bars_Lp_y[pos,3] = min_v[3] + sz * (cont3 - 1)

                        ind_row[starter+num_ele+1] = pos+NAx
                        ind_col[starter+num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]], nodi_red)
                        vals_A[starter+num_ele+1] = -1.0
                        num_ele = num_ele + 1

                        bars_Lp_y[pos,4] = bars_Lp_y[pos,1] + sx
                        bars_Lp_y[pos,5] = bars_Lp_y[pos,2] + sy
                        bars_Lp_y[pos,6] = bars_Lp_y[pos,3] + sz

                        ind_row[starter+num_ele+1] = pos+NAx
                        ind_col[starter+num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont, cont2+1, cont3, Nx, Ny)]], nodi_red)
                        vals_A[starter+num_ele+1] = 1.0
                        num_ele = num_ele + 1

                        liy_mat[pos, 1] = k+1
                        liy_mat[pos, 2] = k+1

                        if cont2 > 1
                           if (matrice[k][cont][cont2-1][cont3]==0)
                               liy_border[pos, 1] = k + 1
                               bars_Lp_y[pos,2] = bars_Lp_y[pos,2] - sy/2.0
                           end
                        else
                            liy_border[pos, 1] = k + 1
                            bars_Lp_y[pos,2] = bars_Lp_y[pos,2] - sy / 2.0
                        end

                        if cont2 + 1 == Ny
                            liy_border[pos, 2] = k + 1
                            bars_Lp_y[pos,5] = bars_Lp_y[pos,5] + sy / 2.0
                        else
                            if cont2 + 2 <= Ny
                                if (matrice[k][cont][cont2+2][cont3]==0)
                                    liy_border[pos, 2] = k + 1
                                    bars_Lp_y[pos,5] = bars_Lp_y[pos,5] + sy / 2.0
                                end
                            end
                        end
                        break
                    end
                end
            end
        end
    end


    liz_mat = zeros(Int8, NAz, 2)
    liz_border = zeros(Int8, NAz, 2)
    bars_Lp_z = zeros(Float64, NAz, 6)

    starter = starter+num_ele
    num_ele = 0

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz - 1)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont][cont2][cont3+1] == 1))

                        pos = mapAz[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)] + 1

                        bars_Lp_z[pos,1] = min_v[1] + sx * (cont - 1)
                        bars_Lp_z[pos,2] = min_v[2] + sy * (cont2 - 1)
                        bars_Lp_z[pos,3] = min_v[3] + sz * (cont3 - 1) + sz/2

                        ind_row[starter+num_ele+1] = pos+NAx+NAy
                        ind_col[starter+num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]], nodi_red)
                        vals_A[starter+num_ele+1] = -1.0
                        num_ele = num_ele + 1

                        bars_Lp_z[pos,4] = bars_Lp_z[pos,1] + sx
                        bars_Lp_z[pos,5] = bars_Lp_z[pos,2] + sy
                        bars_Lp_z[pos,6] = bars_Lp_z[pos,3] + sz

                        ind_row[starter+num_ele+1] = pos+NAx+NAy
                        ind_col[starter+num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont, cont2, cont3+1, Nx, Ny)]], nodi_red)
                        vals_A[starter+num_ele+1] = 1.0
                        num_ele = num_ele + 1

                        liz_mat[pos, 1] = k+1
                        liz_mat[pos, 2] = k+1

                        if cont3 > 1
                           if (matrice[k][cont][cont2][cont3-1]==0)
                               liz_border[pos, 1] = k + 1
                               bars_Lp_z[pos,3] = bars_Lp_z[pos,3] - sz/2.0
                           end
                        else
                            liz_border[pos, 1] = k + 1
                            bars_Lp_z[pos,3] = bars_Lp_z[pos,3] - sz / 2.0
                        end

                        if cont3 + 1 == Nz
                            liz_border[pos, 2] = k + 1
                            bars_Lp_z[pos,6] = bars_Lp_z[pos,6] + sz / 2.0
                        else
                            if cont3 + 2 <= Nz
                                if (matrice[k][cont][cont2][cont3+2]==0)
                                    liz_border[pos, 2] = k + 1
                                    bars_Lp_z[pos,6] = bars_Lp_z[pos,6] + sz / 2.0
                                end
                            end
                        end
                        break
                    end
                end
            end
        end
    end


    return ind_row,ind_col,vals_A,lix_mat,liy_mat,liz_mat,lix_border,liy_border,liz_border,bars_Lp_x,bars_Lp_y,bars_Lp_z
end


#Aggiornata
function compute_diagonals(escalings, materials, sx, sy, sz, lix_mat, liy_mat, liz_mat, lix_border, liy_border, liz_border)
    escaling_R = escalings.R
    escaling_Cd = escalings.Cd
    escaling_Lp = escalings.Lp
    eps0 = 8.854187816997944e-12
    for cont in range(1,length(materials))
        sigmar = materials[cont].sigmar
        epsr = materials[cont].epsr
        if sigmar != 0
            materials[cont].Rx = 0.5 * sx / (sigmar * sy * sz)
            materials[cont].Ry = 0.5 * sy / (sigmar * sx * sz)
            materials[cont].Rz = 0.5 * sz / (sigmar * sy * sx)
            if epsr == 1
                materials[cont].Cx = 0
                materials[cont].Cy = 0
                materials[cont].Cz = 0
            else
                materials[cont].Cx = eps0 * (epsr - 1) * sy * sz / (0.5 * sx)
                materials[cont].Cy = eps0 * (epsr - 1) * sx * sz / (0.5 * sy)
                materials[cont].Cz = eps0 * (epsr - 1) * sy * sx / (0.5 * sz)
            end
        else
            epsr = materials[cont].epsr
            materials[cont].Rx = 0
            materials[cont].Ry = 0
            materials[cont].Rz = 0
            materials[cont].Cx = eps0 * (epsr - 1) * sy * sz / (0.5 * sx)
            materials[cont].Cy = eps0 * (epsr - 1) * sx * sz / (0.5 * sy)
            materials[cont].Cz = eps0 * (epsr - 1) * sy * sx / (0.5 * sz)
        end
    end
    Rx = zeros(size(lix_border,1), 4)
    Ry = zeros(size(liy_border,1), 4)
    Rz = zeros(size(liz_border,1), 4)
    Cx = zeros(size(lix_border,1), 4)
    Cy = zeros(size(liy_border,1), 4)
    Cz = zeros(size(liz_border,1), 4)
    for cont in range(1,length(materials))
        if materials[cont].Rx != 0
            ind_m = findall(x -> x == cont, lix_mat[:, 1])
            Rx[ind_m, 1] .= materials[cont].Rx
            ind_m = findall(x -> x == cont, lix_mat[:, 2])
            Rx[ind_m, 2] .= materials[cont].Rx
            ind_m = findall(x -> x == cont, lix_border[:, 1])
            Rx[ind_m, 3] .= materials[cont].Rx
            ind_m = findall(x -> x == cont, lix_border[:, 2])
            Rx[ind_m, 4] .= +materials[cont].Rx
        end
        if materials[cont].Cx != 0
            ind_m = findall(x -> x == cont, lix_mat[:, 1])
            Cx[ind_m, 1] .= materials[cont].Cx
            ind_m = findall(x -> x == cont, lix_mat[:, 2])
            Cx[ind_m, 2] .= materials[cont].Cx
            ind_m = findall(x -> x == cont, lix_border[:, 1])
            Cx[ind_m, 3] .= materials[cont].Cx
            ind_m = findall(x -> x == cont, lix_border[:, 2])
            Cx[ind_m, 4] .= +materials[cont].Cx
        end
        if materials[cont].Ry != 0
            ind_m = findall(x -> x == cont, liy_mat[:, 1])
            Ry[ind_m, 1] .= materials[cont].Ry
            ind_m = findall(x -> x == cont, liy_mat[:, 2])
            Ry[ind_m, 2] .= materials[cont].Ry
            ind_m = findall(x -> x == cont, liy_border[:, 1])
            Ry[ind_m, 3] .= materials[cont].Ry
            ind_m = findall(x -> x == cont, liy_border[:, 2])
            Ry[ind_m, 4] .= +materials[cont].Ry
        end
        if materials[cont].Cy != 0
            ind_m = findall(x -> x == cont, liy_mat[:, 1])
            Cy[ind_m, 1] .= materials[cont].Cy
            ind_m = findall(x -> x == cont, liy_mat[:, 2])
            Cy[ind_m, 2] .= materials[cont].Cy
            ind_m = findall(x -> x == cont, liy_border[:, 1])
            Cy[ind_m, 3] .= materials[cont].Cy
            ind_m = findall(x -> x == cont, liy_border[:, 2])
            Cy[ind_m, 4] .= +materials[cont].Cy
        end
        if materials[cont].Rz != 0
            ind_m = findall(x -> x == cont, liz_mat[:, 1])
            Rz[ind_m, 1] .= materials[cont].Rz
            ind_m = findall(x -> x == cont, liz_mat[:, 2])
            Rz[ind_m, 2] .= materials[cont].Rz
            ind_m = findall(x -> x == cont, liz_border[:, 1])
            Rz[ind_m, 3] .= materials[cont].Rz
            ind_m = findall(x -> x == cont, liz_border[:, 2])
            Rz[ind_m, 4] .= +materials[cont].Rz
        end
        if materials[cont].Cz != 0
            ind_m = findall(x -> x == cont, liz_mat[:, 1])
            Cz[ind_m, 1] .= materials[cont].Cz
            ind_m = findall(x -> x == cont, liz_mat[:, 2])
            Cz[ind_m, 2] .= materials[cont].Cz
            ind_m = findall(x -> x == cont, liz_border[:, 1])
            Cz[ind_m, 3] .= materials[cont].Cz
            ind_m = findall(x -> x == cont, liz_border[:, 2])
            Cz[ind_m, 4] .= +materials[cont].Cz
        end
    end
    lix_aux = ceil.(Int, lix_border[:, 1] / 100) + ceil.(Int, lix_border[:, 2] / 100)
    liy_aux = ceil.(Int, liy_border[:, 1] / 100) + ceil.(Int, liy_border[:, 2] / 100)
    liz_aux = ceil.(Int, liz_border[:, 1] / 100) + ceil.(Int, liz_border[:, 2] / 100)
    i1x = findall(x -> x == 1, lix_aux)
    i2x = findall(x -> x == 2, lix_aux)
    i1y = findall(x -> x == 1, liy_aux)
    i2y = findall(x -> x == 2, liy_aux)
    i1z = findall(x -> x == 1, liz_aux)
    i2z = findall(x -> x == 2, liz_aux)
    Self_x = Lp_self(sx, sy, sz)
    Self_y = Lp_self(sy, sx, sz)
    Self_z = Lp_self(sz, sy, sx)
    Self_x1 = Lp_self(sx + sx / 2, sy, sz)
    Self_y1 = Lp_self(sy + sy / 2, sx, sz)
    Self_z1 = Lp_self(sz + sz / 2, sy, sx)
    Self_x2 = Lp_self(2 * sx, sy, sz)
    Self_y2 = Lp_self(2 * sy, sx, sz)
    Self_z2 = Lp_self(2 * sz, sy, sx)
    diag_Lp_x = Self_x * ones(size(lix_aux, 1), 1)
    diag_Lp_y = Self_y * ones(size(liy_aux, 1), 1)
    diag_Lp_z = Self_z * ones(size(liz_aux, 1), 1)
    diag_Lp_x[i1x] = Self_x1 * ones(length(i1x), 1)
    diag_Lp_y[i1y] = Self_y1 * ones(length(i1y), 1)
    diag_Lp_z[i1z] = Self_z1 * ones(length(i1z), 1)
    diag_Lp_x[i2x] = Self_x2 * ones(length(i2x), 1)
    diag_Lp_y[i2y] = Self_y2 * ones(length(i2y), 1)
    diag_Lp_z[i2z] = Self_z2 * ones(length(i2z), 1)
    diagonals = Dict()
    diagonals["R"] = escaling_R * [Rx; Ry; Rz]
    diagonals["Cd"] = escaling_Cd * [Cx; Cy; Cz]
    diagonals["Lp"] = escaling_Lp * [diag_Lp_x; diag_Lp_y; diag_Lp_z]
    diagonals["fc_Lp"] = escaling_Lp * ([diag_Lp_x; diag_Lp_y; diag_Lp_z] - [Self_x * ones(size(lix_aux, 1), 1); Self_y * ones(size(liy_aux, 1), 1); Self_z * ones(size(liz_aux, 1), 1)])

    return diagonals
end

#Non sostituita
function create_Gamma_and_center_sup(matrice, Nx,Ny,Nz, map_volumes, min_v, sx, sy, sz, nodi, nodi_red, ext_grids)
    num_grids = size(matrice)[1]
    mapping_surf_1 = zeros(Int64, Nx * Ny * Nz)
    mapping_surf_2 = zeros(Int64, Nx * Ny * Nz)
    mapping_surf_3 = zeros(Int64, Nx * Ny * Nz)
    mapping_surf_4 = zeros(Int64, Nx * Ny * Nz)
    mapping_surf_5 = zeros(Int64, Nx * Ny * Nz)
    mapping_surf_6 = zeros(Int64, Nx * Ny * Nz)

    num_ele_1 = 0
    num_ele_2 = 0
    num_ele_3 = 0
    num_ele_4 = 0
    num_ele_5 = 0
    num_ele_6 = 0

    nnz_surf_max = 6 * Nx * Ny * Nz
    ind_r = ones(Int64, nnz_surf_max)
    ind_c = ones(Int64, nnz_surf_max)

    sup_centers = zeros(Float64, nnz_surf_max,3)
    sup_type = zeros(Int64, nnz_surf_max)


    contat_tot = 1


    cont2=1
    for cont in range(1, stop=Nx)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3]==1) && (ext_grids[k,2,cont,cont2,cont3]==0))
                    num_ele_1 = num_ele_1 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_1[p31] = num_ele_1
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1)  + sx / 2.0
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1)
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                    sup_type[ind_c[contat_tot]] = 1
                    break
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(2, stop=Ny)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3]==1) && (matrice[k][cont][cont2-1][cont3]==0) && (ext_grids[k,2,cont,cont2,cont3]==0))
                        num_ele_1 = num_ele_1 + 1
                        p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                        mapping_surf_1[p31] = num_ele_1
                        contat_tot = contat_tot + 1
                        ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                        ind_c[contat_tot] = num_ele_1
                        sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                        sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1)
                        sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                        sup_type[ind_c[contat_tot]] = 1
                        break
                    end
                end
            end
        end
    end

    cont2 = Ny
    for cont in range(1, stop=Nx)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3] == 1) && (ext_grids[k,1,cont,cont2,cont3] == 0))
                    num_ele_2 = num_ele_2 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_2[p31] = num_ele_2
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1 + num_ele_2
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                    sup_type[ind_c[contat_tot]] = 1
                    break
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny-1)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont][cont2+1][cont3] == 0) && (ext_grids[k,1,cont,cont2,cont3] == 0))
                        check_others = false
                        for k2 in range(1, stop=num_grids)
                            if k!=k2
                                if (matrice[k2][cont][cont2+1][cont3]==1)
                                    check_others = true
                                    break
                                end
                            end
                        end

                        if (check_others == false)
                            num_ele_2 = num_ele_2 + 1
                            p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                            mapping_surf_2[p31] = num_ele_2
                            contat_tot = contat_tot + 1
                            ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                            ind_c[contat_tot] = num_ele_1 + num_ele_2
                            sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                            sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy
                            sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                            sup_type[ind_c[contat_tot]] = 1
                        end
                        break
                    end
                end
            end
        end
    end


    cont = 1
    for cont2 in range(1, stop=Ny)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3] == 1) && (ext_grids[k,4,cont,cont2,cont3] == 0))
                    num_ele_3 = num_ele_3 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_3[p31] = num_ele_3
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1)
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                    sup_type[ind_c[contat_tot]] = 2
                    break
                end
            end
        end
    end

    for cont in range(2, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3]==1) && (matrice[k][cont-1][cont2][cont3]==0) && (ext_grids[k,4,cont,cont2,cont3]==0))
                        num_ele_3 = num_ele_3 + 1
                        p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                        mapping_surf_3[p31] = num_ele_3
                        contat_tot = contat_tot + 1
                        ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                        ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3
                        sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1)
                        sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                        sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                        sup_type[ind_c[contat_tot]] = 2
                        break
                    end
                end
            end
        end
    end

    cont = Nx
    for cont2 in range(1, stop=Ny)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3] == 1) && (ext_grids[k,3,cont,cont2,cont3] == 0))
                    num_ele_4 = num_ele_4 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_4[p31] = num_ele_4
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                    sup_type[ind_c[contat_tot]] = 2
                    break
                end
            end
        end
    end

    for cont in range(1, stop=Nx-1)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont+1][cont2][cont3] == 0) && (ext_grids[k,3,cont,cont2,cont3] == 0))
                        check_others = false
                        for k2 in range(1, stop=num_grids)
                            if k!=k2
                                if (matrice[k2][cont+1][cont2][cont3]==1)
                                    check_others = true
                                    break
                                end
                            end
                        end

                        if (check_others == false)
                            num_ele_4 = num_ele_4 + 1
                            p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                            mapping_surf_4[p31] = num_ele_4
                            contat_tot = contat_tot + 1
                            ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                            ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4
                            sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx
                            sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                            sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                            sup_type[ind_c[contat_tot]] = 2
                        end
                        break
                    end
                end
            end
        end
    end



    cont3 = 1
    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3] == 1) && (ext_grids[k,6,cont,cont2,cont3] == 0))
                    num_ele_5 = num_ele_5 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_5[p31] = num_ele_5
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4 + num_ele_5
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1)
                    sup_type[ind_c[contat_tot]] = 3
                    break
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(2, stop=Nz)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3]==1) && (matrice[k][cont][cont2][cont3-1]==0) && (ext_grids[k,6,cont,cont2,cont3]==0))
                        num_ele_5 = num_ele_5 + 1
                        p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                        mapping_surf_5[p31] = num_ele_5
                        contat_tot = contat_tot + 1
                        ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                        ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4 + num_ele_5
                        sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                        sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                        sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1)
                        sup_type[ind_c[contat_tot]] = 3
                        break
                    end
                end
            end
        end
    end

    cont3 = Nz
    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3] == 1) && (ext_grids[k,5,cont,cont2,cont3] == 0))
                    num_ele_6 = num_ele_6 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_6[p31] = num_ele_6
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4 + num_ele_5 + num_ele_6
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz
                    sup_type[ind_c[contat_tot]] = 3
                    break
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz-1)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont][cont2][cont3+1] == 0) && (ext_grids[k,5,cont,cont2,cont3] == 0))
                        check_others = false
                        for k2 in range(1, stop=num_grids)
                            if k!=k2
                                if (matrice[k2][cont][cont2][cont3+1]==1)
                                    check_others = true
                                    break
                                end
                            end
                        end

                        if (check_others==false)
                            num_ele_6 = num_ele_6 + 1
                            p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                            mapping_surf_6[p31] = num_ele_6
                            contat_tot = contat_tot + 1
                            ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                            ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4 + num_ele_5 + num_ele_6
                            sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                            sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                            sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz
                            sup_type[ind_c[contat_tot]] = 3
                        end
                        break
                    end
                end
            end
        end
    end

    nnz_surf = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4 + num_ele_5 + num_ele_6
    vals = ones(Float64, contat_tot)
    ind_r = ind_r[1:contat_tot]
    ind_c = ind_c[1:contat_tot]
    sup_centers = sup_centers[1:nnz_surf,:]
    sup_type = sup_type[1:nnz_surf]


    return vals,ind_r,ind_c,sup_centers,sup_type
end


function generate_interconnection_matrices_and_centers(escalings,size_x,size_y,size_z,grid_matrix,num_cel_x,num_cel_y,num_cel_z,materials,port_matrix,lumped_el_matrix,minimum_vertex)

    dominant_list=1;

    numTotVox = num_cel_x*num_cel_y*num_cel_z
    
    println("Total Number of Voxels (including air):", numTotVox)
    n_grids = length(materials)
    @assert size(grid_matrix)[1]==n_grids
    
    num_tot_full_vox = 0

    
    for i in range(1, stop=n_grids)
        for j in range(1, stop=num_cel_x)
            for k in range(1, stop=num_cel_y)
                num_tot_full_vox = num_tot_full_vox + count(i-> i==1,grid_matrix[i][j][k])
            end
        end
    end

    # TODO: VERIFICARE SE OCCORRE
    cont=1
    while cont<=length(materials)
        materials[cont].epsr = materials[cont].permittivity +1im * materials[cont].permittivity * materials[cont].tangent_delta_permittivity
        cont+=1
    end

    @assert num_cel_x isa Int
    @assert num_cel_y isa Int
    @assert num_cel_z isa Int

    n_boxes,mapping_vols,volume_centers,volumes_materials=create_volumes_mapping_and_centers(grid_matrix,num_cel_x,num_cel_y,num_cel_z,num_tot_full_vox,size_x,size_y,size_z,minimum_vertex)

    externals_grids = create_external_grids(grid_matrix,num_cel_x,num_cel_y,num_cel_z)
    
    nodes_red, nodes = create_nodes_ref(grid_matrix,num_cel_x,num_cel_y,num_cel_z,num_tot_full_vox,externals_grids,mapping_vols, dominant_list)

    port_matrix.port_voxels, port_matrix.port_nodes = find_nodes_port(volume_centers, port_matrix.port_start, port_matrix.port_end, nodes, nodes_red)
        

    lumped_el_matrix.le_voxels, lumped_el_matrix.le_nodes = find_nodes_port(volume_centers,lumped_el_matrix.le_start,lumped_el_matrix.le_end, nodes, nodes_red)

    vals,ind_r,ind_c,sup_centers,sup_type=create_Gamma_and_center_sup(grid_matrix,num_cel_x,num_cel_y,num_cel_z,mapping_vols, minimum_vertex, size_x,size_y,size_z, nodes, nodes_red, externals_grids)

    Gamma = sparse(ind_r, ind_c, vals)

    println("Number of Surfaces (without air):", size(Gamma)[2])

    map_for_Ax, n_for_Ax = create_mapping_Ax(grid_matrix,mapping_vols, nodes, nodes_red)
    map_for_Ay, n_for_Ay = create_mapping_Ay(grid_matrix,mapping_vols, nodes, nodes_red)
    map_for_Az, n_for_Az = create_mapping_Az(grid_matrix,mapping_vols, nodes, nodes_red)

    ind_row,ind_col,vals_A, LiX, LiY, LiZ, LiX_bord,LiY_bord,LiZ_bord,bars_Lp_x,bars_Lp_y,bars_Lp_z = create_A_mats_volInd(grid_matrix,num_cel_x,num_cel_y,num_cel_z,mapping_vols,
                                       map_for_Ax, n_for_Ax, map_for_Ay, n_for_Ay, map_for_Az, n_for_Az,
                                       size_x,size_y,size_z,minimum_vertex,nodes,nodes_red)

    

    A = sparse(ind_row, ind_col, vals_A)

    println("Edges without air:", n_for_Ax+n_for_Ay+n_for_Az)
    println("Nodes without air:", size(Gamma)[1])

    diagonals=compute_diagonals(escalings, materials, size_x,size_y,size_z, LiX, LiY, LiZ,LiX_bord, LiY_bord, LiZ_bord)
    

    return A, Gamma, port_matrix, lumped_el_matrix, sup_centers, sup_type, bars_Lp_x, bars_Lp_y, bars_Lp_z, diagonals["R"], diagonals["Cd"]
end