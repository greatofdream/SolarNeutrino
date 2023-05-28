.PHONY: all
# Data
reactions:=N13 O15 F17 pp B8 hep
BP2004: Data/BP2004/gs98/model.dat Data/BP2004/gs98/flux.dat
BSBGS98: Data/BP2005/gs98/model.dat Data/BP2005/gs98/flux.dat Data/BP2005/gs98/preview.h5
BSBAGS05: Data/BP2005/ags05/model.dat Data/BP2005/ags05/flux.dat
SPECTRA: Data/SPECTRA/Be7.dat Data/SPECTRA/B8.dat Data/SPECTRA/N13.dat Data/SPECTRA/O15.dat Data/SPECTRA/F17.dat Data/SPECTRA/pp.dat Data/SPECTRA/hep.dat Data/SPECTRA/preview.h5
CROSSSECTION:
# Data
Data/SPECTRA/Be7.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/7Belineshape/7belineshape.dat -O $@
Data/SPECTRA/B8.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/B8spectrum/b8spectrum.txt -O $@
Data/SPECTRA/N13.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/n13.dat -O $@
Data/SPECTRA/O15.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/o15.dat -O $@
Data/SPECTRA/F17.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/CNOspectra/f17.dat -O $@
Data/SPECTRA/pp.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/PPenergyspectrum/ppenergytab -O $@
Data/SPECTRA/hep.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/Hepspectrum/hepspectrum.dat -O $@
Data/SPECTRA/T_table.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/Momentsspectra/T_table.dat -O $@
Data/SPECTRA/T2_table.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/Momentsspectra/T2_table.dat -O $@
Data/SPECTRA/preview.h5: $(reactions:%=Data/SPECTRA/%.dat)
	python3 SpectraPreview.py -i $^ --reactions $(reactions) -o $@ --models Data/BP2005/gs98/preview.h5
Data/BP2004/gs98/model.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BP2004/bp2004stdmodel.dat -O $@
Data/BP2004/gs98/flux.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BP2004/bp2004flux.dat -O $@

Data/BP2005/gs98/model.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs05op.dat -O $@
Data/BP2005/gs98/flux.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs2005opflux.dat -O $@

Data/BP2005/ags05/model.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs05_agsop.dat -O $@
Data/BP2005/ags05/flux.dat:
	mkdir -p $(dir $@)
	wget http://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs2005agsopflux.dat -O $@
# Data Preview
Data/BP%/preview.h5: Data/BP%/model.dat  Data/BP%/flux.dat
	python3 SSMPreview.py -i $^ -o $@ --ssm BSB > $@.log 2>&1
	python3 SSMPreview.py -i $^ -o $@ --ssm BSB --plot
Data/Mesa%/preview.h5: Data/Mesa%/model.dat
	python3 SSMPreview.py -i $^ -o $@ --ssm Mesa > $@.log 2>&1
	python3 SSMPreview.py -i $^ -o $@ --ssm Mesa --plot 
Data/B16%/preview.h5: Data/B16%/model.dat
	python3 SSMPreview.py -i $^ -o $@ --ssm B16 > $@.log 2>&1
	python3 SSMPreview.py -i $^ -o $@ --ssm B16 --plot 
# Predict
Predict/%/fluxSolar.h5: Data/%/preview.h5
	python3 evolutionSolar.py -i $^ -o $@
Predict/%/fluxVaccum.h5: Predict/%/fluxSolar.h5
	python3 evolutionVaccum.py -i $^ -o $@
Predict/%/fluxEarth.h5: Predict/%/fluxVaccum.h5
	python3 evolutionEarth.py -i $@ -o $@
Predict/%.pdf: Predict/%.h5
	python3 evolutionPlot.py -i $^ -o $@
