module Sz
    function r(matrix)
        # returns size of matrix as (rows, columns)
        mat_pair = size(matrix)[1]
        return mat_pair
    end  

    function c(matrix)
        mat_pair = size(matrix)[2]
        return mat_pair
    end

    function z(matrix)
        z = (r(matrix) == 0)
    end
end
