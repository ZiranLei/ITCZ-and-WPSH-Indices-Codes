; ==============================================================================
; Script: calc_ITCZ_lat.ncl
; 
; Description: This script calculates the mean latitudinal position of the ITCZ 
;              based on ERA5 500-hPa vertical velocity.
; 
; Input Data:  Vertical_velocity.nc (Variable: w)
; 
; Output:      ITCZ_lat_output.csv
; 
; Requirements: NCAR Command Language (NCL) Version 6.6.2 or later.
; ==============================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
    ; ---Read Data ---
    f = addfile("Vertical_velocity.nc", "r")
    omega = tofloat(f->w)
    valid_time = f->valid_time
    latitude = f->latitude
    longitude = f->longitude
    
    ; ---Define Regions and Thresholds ---
    lon_min = 90.0
    lon_max = 130.0
    lat_min = -20.0
    lat_max = 20.0
    
    lat_idx_int = tointeger(ind(latitude.ge.lat_min .and. latitude.le.lat_max))
    lon_idx_int = tointeger(ind(longitude.ge.lon_min .and. longitude.le.lon_max))
    
    lat_sub = latitude(lat_idx_int)
    lon_sub = longitude(lon_idx_int)
    
    ntime     = dimsizes(valid_time)
    nlat_sub  = dimsizes(lat_sub)
    nlon_sub  = dimsizes(lon_sub)
    
    ; ---Time Formatting ---
    time_float = todouble(valid_time)
    time_float@units = valid_time@units
    time_float@calendar = valid_time@calendar
    time_cal = cd_calendar(time_float, 0)
    
    years  = toint(time_cal(:,0))
    months = toint(time_cal(:,1))
    days   = toint(time_cal(:,2))
    
    dates = new(ntime, "string")
    do t = 0, ntime-1
        dates(t) = sprinti("%04d", years(t)) + "-" + sprinti("%02d", months(t)) + "-" + sprinti("%02d", days(t))
    end do
    
    ; ---Initialize Arrays ---
    ITCZ_lat = new(ntime, "float")
    ITCZ_lat@_FillValue = -999.0
    
    print("Processing " + ntime + " time steps for ITCZ latitude...")
    
    ; ---Main Calculation Loop ---
    do t = 0, ntime-1
        if (mod(t, 100) .eq. 0) then
            print("   Processing time step: " + t + " / " + ntime)
        end if
        
        weighted_sum = 0.0d0
        total_weight = 0.0d0
        
        do ii = 0, nlat_sub-1
            ilat = lat_idx_int(ii)
            current_lat = lat_sub(ii)
            
            do jj = 0, nlon_sub-1
                ilon = lon_idx_int(jj)
                omega_val = omega(t, 0, ilat, ilon)  ; time, plev, lat, lon
                
                ; Filter missing values and descending motion (keep only ω < 0)
                if (ismissing(omega_val) .or. omega_val .ge. 0.0) then
                    continue
                end if
                
                ; Calculate weight (|ω|)
                weight = abs(omega_val)
                
                ; Accumulate weighted sum and total weight
                weighted_sum = weighted_sum + (current_lat * weight)
                total_weight = total_weight + weight
            end do
        end do
        
        ; Calculate the weighted mean latitude
        if (total_weight .gt. 0.0) then
            ITCZ_lat(t) = tofloat(weighted_sum / total_weight)
        else
            ITCZ_lat(t) = ITCZ_lat@_FillValue
        end if
    end do
    
    ; --- Output CSV Files ---
    out_csv = "ITCZ_lat_output.csv"
    if (fileexists(out_csv)) then
        system("rm -f " + out_csv)
    end if
    
    ; Write header
    header = "Date,ITCZ_Latitude(degrees_north)"
    write_table(out_csv, "w", [/header/], "%s")
    
    ; Write data lines
    do t = 0, ntime-1
        if (.not. ismissing(ITCZ_lat(t))) then
            lat_str = sprintf("%.3f", ITCZ_lat(t))
        else
            lat_str = "NA"
        end if
        
        line = dates(t) + "," + lat_str
        write_table(out_csv, "a", [/line/], "%s")
    end do
    
    print("Processing complete. " + out_csv + " generated.")
end