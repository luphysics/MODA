function x = csv_to_mvar( s )
%CSV_TO_MVAR converts comma separated string to a vector
x = eval( [ '[', s, ']' ] );
end

