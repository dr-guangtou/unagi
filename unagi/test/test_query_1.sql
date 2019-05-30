SELECT 
        object_id 
      , ra 
      , dec 
      , imag_kron 
      , imag_kron_err 
      , ymag_kron 
      , ymag_kron_err 
      , imag_kron - ymag_kron AS i_y 
    FROM 
        pdr1_udeep.forced 
    WHERE 
        boxsearch(coord, 34.0, 36.0, -5.0, -4.5) 
        /* is equivalent to 
             ra BETWEEN 34.0 AND 36.0 
             AND dec BETWEEN -5.0 AND -4.5 
           but boxSearch() is much faster 
        */ 
        AND imag_kron < 25.5 
    LIMIT 10
;