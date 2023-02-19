

band='u'

cluster_name=cl_m4.1e14_z0.059

for dir in $cluster_name/r{0..29}; do
    yaml_file=$dir/forecast_"$band"_"$cluster_name".yaml
    gs_file=$dir/forecast_"$band"_gs_config.yaml
  if [ -f $yaml_file ]; then
      sed -i -e 's/#- galsim/- galsim/g; s/- medsmaker/#- medsmaker/g; s/- metacal/#- metacal/g; s/- shear_profile/#- shear_profile/g ' \
	  -e 's/overwrite: false/overwrite: true/g' $yaml_file
      sed -i 's/nexp: 36/nexp: 1/g' $gs_file
  fi
done

cluster_name=cl_m4.1e14_z0.3

for dir in $cluster_name/r{0..29}; do
    yaml_file=$dir/forecast_"$band"_"$cluster_name".yaml
    gs_file=$dir/forecast_"$band"_gs_config.yaml
  if [ -f $yaml_file ]; then
      sed -i -e 's/#- galsim/- galsim/g; s/- medsmaker/#- medsmaker/g; s/- metacal/#- metacal/g; s/- shear_profile/#- shear_profile/g ' \
	  -e 's/overwrite: false/overwrite: true/g' $yaml_file
      sed -i 's/nexp: 36/nexp: 1/g' $gs_file
  fi
done

cluster_name=cl_m4.1e14_z0.45

for dir in $cluster_name/r{0..29}; do
    yaml_file=$dir/forecast_"$band"_"$cluster_name".yaml
    gs_file=$dir/forecast_"$band"_gs_config.yaml
  if [ -f $yaml_file ]; then
      sed -i -e 's/#- galsim/- galsim/g; s/- medsmaker/#- medsmaker/g; s/- metacal/#- metacal/g; s/- shear_profile/#- shear_profile/g ' \
	  -e 's/overwrite: false/overwrite: true/g' $yaml_file
      sed -i 's/nexp: 36/nexp: 1/g' $gs_file
  fi
done



#########

cluster_name=cl_m4.1e14_z0.3

for dir in $cluster_name/r{0..29}; do
    yaml_file=$dir/forecast_"$band"_"$cluster_name".yaml
  if [ -f $yaml_file ]; then
      sed -i -e 's/#- galsim/- galsim/g' $yaml_file
  fi
done
