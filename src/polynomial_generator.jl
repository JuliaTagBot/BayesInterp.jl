struct poly_exponents{N,T}
    data::Matrix{T}
end





@generated function max_degree(::Val{N}, D) where N
    quote
        out = Matrix{Int}(undef, $N, binomial($N+D,$N))
        r_0 = D
        ind = 0
        Base.Cartesian.@nloops $N i j -> begin
            0:r_{$N-j}            
        end j -> begin
            r_{$(N+1)-j} = r_{$N-j} - i_j
        end begin
            ind += 1
            Base.Cartesian.@nexprs $N j -> ( out[j,ind] = i_j )
        end
        poly_exponents{$N,Int}(out)
    end
end

function norm_degree_quote(N::Int)
    Np1 = N+1
    @static if VERSION > v"0.6.9" # should be PR that changed reinterpret
        outgen = quote
            out = Matrix{Int}(undef, $N, length(outc))
            @inbounds for i = eachindex(outc)
                out[:,i] .= outc[i]
            end
        end
    else
        outgen = :( out = reshape(reinterpret( Int, outc ), $N, length(outc)) )
    end
    quote
        Li = inv(L)
        j_0 = D
        s_0 = 0
        DL = D^L
        outc = Vector{NTuple{$N,Int}}(undef, 0)
        sizehint!(outc, round(Int,3*($N*D^L)^Li)) 
        @nloops $N i p -> begin
            0:j_{$N-p}
        end p -> begin
            iL_p = i_p^L
            s_{$Np1-p} = s_{$N-p} + i_p^L
            j_{$Np1-p} = floor(Int, (DL - s_{$Np1-p})^Li)
        end begin
            push!(outc, @ntuple $N i )
        end
        $outgen
        poly_exponents{$N,Int}(out)
    end
end
@generated function norm_degree(::Val{N}, D, L) where N
    norm_degree_quote(N)
end




@generated function max_degree_vec(::Val{N}, D) where N
    quote
        out = Vector{NTuple{$N,Int}}(undef, binomial($N+D,$N))
        r_0 = D
        ind = 0
        Base.Cartesian.@nloops $N i j -> begin
            0:r_{$N-j}            
        end j -> begin
            r_{$(N+1)-j} = r_{$N-j} - i_j
        end begin
            ind += 1
            out[ind] = Base.Cartesian.@ntuple $N i
        end
        out
    end
end

@generated function norm_degree_vec(::Val{N}, D, L) where N
    Np1 = N+1
    quote
        Li = inv(L)
        j_0 = D
        DL = D^L
        s_0 = 0
        out = Vector{NTuple{$N,Int}}(undef, 0)
        @nloops $N i p -> begin
            0:j_{$N-p}
        end p -> begin
            s_{$Np1-p} = s_{$N-p} + i_p^L
            j_{$Np1-p} = floor(Int, (DL - s_{$Np1-p})^Li)
        end begin
            push!(out, @ntuple $N i )
        end
        out
    end
end



