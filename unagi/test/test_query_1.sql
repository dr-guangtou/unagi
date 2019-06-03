SELECT
	object_id
      , ra
      , dec
      , i_kronflux_mag
      , i_kronflux_magsigma
      , y_kronflux_mag
      , y_kronflux_magsigma
      , i_kronflux_mag - y_kronflux_mag AS i_y
    FROM
	pdr2_dud.forced
	JOIN pdr2_dud.forced2 USING (object_id)
    WHERE
	  boxSearch(coord, 34.0, 36.0, -5.0, -4.5)
          /* is equivalent to
                 ra  BETWEEN 34.0 AND 36.0
             AND dec BETWEEN -5.0 AND -4.5
             but boxSearch() is much faster
          */
	AND i_kronflux_mag < 25.5
    LIMIT 10
;