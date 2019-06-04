-- Get "clean" objects using several flags
-- Works for dr2.s18a and pdr2

SELECT
        object_id
        ,f1.ra
        ,f1.dec
        ,f1.g_cmodel_flux
        ,f1.g_cmodel_fluxsigma
        ,f1.r_cmodel_flux
        ,f1.r_cmodel_fluxsigma
        ,f1.i_cmodel_flux
        ,f1.i_cmodel_fluxsigma
        ,f1.z_cmodel_flux
        ,f1.z_cmodel_fluxsigma
        ,f1.y_cmodel_flux
        ,f1.y_cmodel_fluxsigma
        ,f2.g_psfflux_flux
        ,f2.g_psfflux_fluxsigma
        ,f2.r_psfflux_flux
        ,f2.r_psfflux_fluxsigma
        ,f2.i_psfflux_flux
        ,f2.i_psfflux_fluxsigma
        ,f2.z_psfflux_flux
        ,f2.z_psfflux_fluxsigma
        ,f2.y_psfflux_flux
        ,f2.y_psfflux_fluxsigma
        ,f1.i_extendedness_value
        ,f1.a_g
        ,f1.a_r
        ,f1.a_i
        ,f1.a_z
        ,f1.a_y
    FROM
        {RERUN_USED}.forced  AS f1
      LEFT JOIN
        {RERUN_USED}.forced2 AS f2 USING (object_id)
      LEFT JOIN
        {RERUN_USED}.meas   AS m1 USING (object_id)
      LEFT JOIN
        {RERUN_USED}.meas2  AS m2 USING (object_id)
    WHERE
        f1.isprimary
        AND NOT m2.g_sdsscentroid_flag
        AND NOT m2.r_sdsscentroid_flag
        AND NOT m2.i_sdsscentroid_flag
        AND NOT m2.z_sdsscentroid_flag
        AND NOT m2.y_sdsscentroid_flag
        AND NOT f1.g_pixelflags_edge
        AND NOT f1.r_pixelflags_edge
        AND NOT f1.i_pixelflags_edge
        AND NOT f1.z_pixelflags_edge
        AND NOT f1.y_pixelflags_edge
        AND NOT f1.g_pixelflags_interpolatedcenter
        AND NOT f1.r_pixelflags_interpolatedcenter
        AND NOT f1.i_pixelflags_interpolatedcenter
        AND NOT f1.z_pixelflags_interpolatedcenter
        AND NOT f1.y_pixelflags_interpolatedcenter
        AND NOT f1.g_pixelflags_saturatedcenter
        AND NOT f1.r_pixelflags_saturatedcenter
        AND NOT f1.i_pixelflags_saturatedcenter
        AND NOT f1.z_pixelflags_saturatedcenter
        AND NOT f1.y_pixelflags_saturatedcenter
        AND NOT f1.g_pixelflags_crcenter
        AND NOT f1.r_pixelflags_crcenter
        AND NOT f1.i_pixelflags_crcenter
        AND NOT f1.z_pixelflags_crcenter
        AND NOT f1.y_pixelflags_crcenter
        AND NOT f1.g_pixelflags_bad
        AND NOT f1.r_pixelflags_bad
        AND NOT f1.i_pixelflags_bad
        AND NOT f1.z_pixelflags_bad
        AND NOT f1.y_pixelflags_bad
        AND NOT f1.g_cmodel_flag
        AND NOT f1.r_cmodel_flag
        AND NOT f1.i_cmodel_flag
        AND NOT f1.z_cmodel_flag
        AND NOT f1.y_cmodel_flag
;