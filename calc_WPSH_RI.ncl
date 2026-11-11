; ==============================================================================
; Script: calc_WPSH_RI.ncl
; 
; Description: This script calculates the Western Pacific Subtropical High 
;              Ridge Index (WPSHRI) based on ERA5 500-hPa Geopotential
;              and zonal wind fields.
; 
; Input Data: 1. Geopotential.nc (Variable: z)
;                      2. U_component_wind.nc   (Variable: u)
;     
; Output:     WPSH_Ridge_Line_Index.csv
; 
; Requirements: NCAR Command Language (NCL) Version 6.6.2 or later.
; ==============================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
    ; --- Read Data ---
    f_z = addfile("Geopotential.nc", "r")
    f_u = addfile("U_component_wind.nc", "r")
    
    z = tofloat(f_z->z)
    u = tofloat(f_u->u)
    
    lat = f_z->latitude
    lon = f_z->longitude
    
    z_ndims = dimsizes(dimsizes(z))
    z_gpm = z / tofloat(9.80665)
    
    ; --- Define Regions and Thresholds ---
    lat_min = 10.0
    lat_max = 60.0
    lon_min = 110.0
    lon_max = 150.0
    z_threshold = 5880.0
    
    lat_idx_int = tointeger(ind(lat.ge.lat_min .and. lat.le.lat_max))
    lon_idx_int = tointeger(ind(lon.ge.lon_min .and. lon.le.lon_max))
    
    ; --- Extract Regional Data ---
    if (z_ndims .eq. 3) then
        z_region = z_gpm(:, lat_idx_int, lon_idx_int)
        u_region = u(:, lat_idx_int, lon_idx_int)
    else if (z_ndims .eq. 4) then
        z_region = z_gpm(:, 0, lat_idx_int, lon_idx_int)
        u_region = u(:, 0, lat_idx_int, lon_idx_int)
    end if
    end if
    
    lat_region = lat(lat_idx_int)
    lon_region = lon(lon_idx_int)
    
    ntime = tointeger(dimsizes(z_region(:,0,0)))
    nlat  = tointeger(dimsizes(lat_region))
    nlon  = tointeger(dimsizes(lon_region))
    
    ; ---  Initialize Arrays ---
    ridge_lat = new((/ntime, nlon/), "float")
    ridge_lat@_FillValue = -999.0
    
    valid_ridge_count = new(ntime, "integer")
    
    avg_ridge_lat = new(ntime, "float")
    avg_ridge_lat@_FillValue = -999.0
    
    print("Processing " + ntime + " time steps for WPSHRI...")
    
    ; --- Main Calculation Loop ---
    do t = 0, ntime-1
        if (mod(t, 100) .eq. 0) then
            print("   Processing time step: " + t + " / " + ntime)
        end if
        
        valid_ridge_count(t) = 0
        sum_ridge_lat = 0.0
        
        ; Step 1: Find WPSH Ridge Line
        do i = 0, nlon-1
            u_profile = u_region(t, :, i)
            z_profile = z_region(t, :, i)
            found_ridge = False
            
            do j = 1, nlat-1
                ; 过滤缺测值
                if (ismissing(u_profile(j-1)) .or. ismissing(u_profile(j)) .or. \
                    ismissing(z_profile(j-1))   .or. ismissing(z_profile(j))) then
                    continue
                end if
                
                ; 过滤未跨越 0 的风向，或跨度太小的无效点
                if ((u_profile(j-1) * u_profile(j) .gt. 0) .or. \
                    (abs(u_profile(j) - u_profile(j-1)) .le. 0.0001)) then
                    continue
                end if
                
                ; 核心计算：线性插值
                weight = -u_profile(j-1) / (u_profile(j) - u_profile(j-1))
                interp_lat = tofloat(lat_region(j-1) * (1.0 - weight) + lat_region(j) * weight)
                interp_z   = tofloat(z_profile(j-1)  * (1.0 - weight) + z_profile(j)  * weight)
                
                ; 检查插值点的位势高度是否 >= 5880
                if (interp_z .ge. z_threshold) then
                    ridge_lat(t, i) = interp_lat
                    sum_ridge_lat = sum_ridge_lat + interp_lat
                    valid_ridge_count(t) = valid_ridge_count(t) + 1
                    found_ridge = True
                    break
                end if
            end do
            
            if (.not.found_ridge) then
                ridge_lat(t, i) = ridge_lat@_FillValue
            end if
        end do
        
        ; Step 2: Calculate Averages
        if (valid_ridge_count(t) .gt. 0) then
            avg_ridge_lat(t) = sum_ridge_lat / tofloat(valid_ridge_count(t))
        else
            avg_ridge_lat(t) = avg_ridge_lat@_FillValue
        end if
    end do
    
    ; --- Time Formatting ---
    time_str = new(ntime, "string")
    if (isfilevar(f_z, "time")) then
        time_var = f_z->time
        if (isatt(time_var, "units")) then
            time_str = cd_string(time_var, "%Y-%m-%d")
        else
            do t = 0, ntime-1
                time_str(t) = "Step_" + tostring(t+1)
            end do
        end if
    else
        do t = 0, ntime-1
            time_str(t) = "Step_" + tostring(t+1)
        end do
    end if
    
    ; --- Output CSV Files ---
    ridge_filename = "WPSH_Ridge_Line_Index.csv"
    if (fileexists(ridge_filename)) then
        system("rm -f " + ridge_filename)
    end if
    
    ridge_header = "Time,Avg_Ridge_Latitude(N),Valid_Ridge_Count"
    do i = 0, nlon-1
        ridge_header = ridge_header + ",Lon_" + sprintf("%.1f", lon_region(i)) + "_Ridge"
    end do
    write_table(ridge_filename, "w", [/ridge_header/], "%s")
    
    do t = 0, ntime-1
        avg_lat_str = "NA"
        if (.not.ismissing(avg_ridge_lat(t))) then
            avg_lat_str = sprintf("%.3f", avg_ridge_lat(t))
        end if
        
        ridge_line = time_str(t) + "," + avg_lat_str + "," + tostring(valid_ridge_count(t))
        do i = 0, nlon-1
            if (.not.ismissing(ridge_lat(t, i))) then
                ridge_line = ridge_line + "," + sprintf("%.3f", ridge_lat(t, i))
            else
                ridge_line = ridge_line + ",NA"
            end if
        end do
        write_table(ridge_filename, "a", [/ridge_line/], "%s")
    end do
    
    print("Processing complete. WPSH_Ridge_Line_Index.csv generated.")
end
