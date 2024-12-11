.PHONY: all
# data
reactions:=N13 O15 F17 pp B8 hep Be7
# BP2004
BP2004: data/BP2004/gs98/model.dat data/BP2004/gs98/flux.dat
# BP2005
BSBGS98: data/BP2005/gs98/model.dat data/BP2005/gs98/flux.dat data/BP2005/gs98/preview.h5
BSBAGS05: data/BP2005/ags05/model.dat data/BP2005/ags05/flux.dat
# B16
B16: data/B16/gs98/model.dat data/B16/gs98/flux.dat data/B16/agss09/model.dat data/B16/agss09/flux.dat data/B16/spectra.dat
# X25
X25GS98:
X25AGSS09:
X25MB22:

# spectra shpae
SPECTRA: data/SPECTRA/Be7.dat data/SPECTRA/B8.dat data/SPECTRA/N13.dat data/SPECTRA/O15.dat data/SPECTRA/F17.dat data/SPECTRA/pp.dat data/SPECTRA/hep.dat data/SPECTRA/preview.h5
CROSSSECTION:
OSCPARAS: data/OSCPARAS/v52.release-SKoff-NO.txt data/OSCPARAS/v52.release-SKoff-IO.txt
# data
data/OSCPARAS/%:
	mkdir -p $(dir $@)
	wget http://www.nu-fit.org/sites/default/files/$*.xz -O $@
	unxz -d $@.xz
data/SPECTRA/Be7.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/7Belineshape/7belineshape.dat -O $@
data/SPECTRA/B8.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/B8spectrum/b8spectrum.txt -O $@
data/SPECTRA/N13.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/n13.dat -O $@
data/SPECTRA/O15.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/o15.dat -O $@
data/SPECTRA/F17.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/f17.dat -O $@
data/SPECTRA/pp.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/PPenergyspectrum/ppenergytab -O $@
data/SPECTRA/hep.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/Hepspectrum/hepspectrum.dat -O $@
data/SPECTRA/T_table.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/Momentsspectra/T_table.dat -O $@
data/SPECTRA/T2_table.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/Momentsspectra/T2_table.dat -O $@
data/SPECTRA/preview.h5: $(reactions:%=data/SPECTRA/%.dat)
	python3 SpectraPreview.py -i $^ pep.dat --reactions $(reactions) pep -o $@ --models data/BP2005/gs98/preview.h5

data/BP2004/gs98/model.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BP2004/bp2004stdmodel.dat -O $@
data/BP2004/gs98/flux.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BP2004/bp2004flux.dat -O $@

data/BP2005/gs98/model.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs05op.dat -O $@
data/BP2005/gs98/flux.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs2005opflux.dat -O $@

data/BP2005/ags05/model.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs05_agsop.dat -O $@
data/BP2005/ags05/flux.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs2005agsopflux.dat -O $@

data/B16/gs98/model.dat:
	mkdir -p $(dir $@)
	wget https://aliga.ice.csic.es/personal/aldos/Solar_Data_files/struct_b16_gs98.dat -O $@
data/B16/gs98/flux.dat:
	mkdir -p $(dir $@)
	wget https://aliga.ice.csic.es/personal/aldos/Solar_Data_files/nudistr_b16_gs98.dat -O $@

data/B16/agss09/model.dat:
	mkdir -p $(dir $@)
	wget https://aliga.ice.csic.es/personal/aldos/Solar_Data_files/struct_b16_agss09.dat -O $@
data/B16/agss09/flux.dat:
	mkdir -p $(dir $@)
	wget https://aliga.ice.csic.es/personal/aldos/Solar_Data_files/nudistr_b16_agss09.dat -O $@

data/B16/flux.dat:
	mkdir -p $(dir $@)
	wget https://aliga.ice.csic.es/personal/aldos/Solar_Data_files/fluxes_b16.dat -O $@

# data Preview
data/%/preview.h5: data/%/model.dat  data/%/flux.dat
	python3 SSMPreview.py -i $^ -o $@ --ssm $(*) > $@.log 2>&1
	python3 SSMPreview.py -i $@ -o $@.pdf --ssm $(*) --plot
data/Mesa%/preview.h5: data/Mesa%/model.dat
	python3 SSMPreview.py -i $^ -o $@ --ssm Mesa > $@.log 2>&1
	python3 SSMPreview.py -i $^ -o $@ --ssm Mesa --plot 
# Predict
data/%/fluxSolar.h5: data/%/preview.h5 $(reactions:%=data/SPECTRA/%.dat)
	python3 NeutrinoEvolution.py -i $< -o $@ --spectra $(reactions:%=data/SPECTRA/%.dat) --reactions $(reactions)
	python3 NeutrinoEvolution.py -i $@ -o $@.pdf --reactions $(reactions) pep --plot
data/%/fluxEarth.h5: Predict/%/fluxVaccum.h5
	python3 evolutionEarth.py -i $@ -o $@
