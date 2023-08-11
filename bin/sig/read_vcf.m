function SVstruct_dd = read_vcf( SVstruct,files )

warning off

% generate SV table for signficance analysis

%DIR=dir(strcat(SV_DIR,'/*.sv.vcf'));
%files={DIR.name};
num_vcf=length(files);

%SVstruct=table();
gc=height(SVstruct)+1;
for c1=1:num_vcf,
    %disp(['...parsing file ' num2str(c1) ' out of ' num2str(num_vcf)]);
    %VCF=importdata(strcat(SV_DIR,'/',files{c1}),'\t');
    VCF=importdata(files{c1},'\t');
    cline=1;

    while isempty(regexp(VCF{cline,1},'##SAMPLE=<ID='))
        cline=cline+1;
    end
    %sampleID=regexp(VCF{cline},'(?<=\,SampleName=)(?:\w|-)*','match');
    %sampleID=regexp(VCF{cline},'(?:\w|-|\.)*(?=.bam)','match','once');
sampleID=regexp(files{c1},'(?:\w|-|\.)*(?=.snowman)','match','once')

    while isempty(regexp(VCF{cline,1},'#CHROM'))
        cline=cline+1;
    end

    for c2=cline+1:length(VCF)

	line=VCF{c2};
        SVstruct.sid{gc,1}=sampleID;
        SVstruct.donor_unique_id{gc,1}=sampleID;
        tokens=regexp(line,'(\w+)\t(\d+)\t(\w+)','tokens','once');
        SVstruct.seqnames(gc,1)=chr_str2num(tokens(1));
        SVstruct.start(gc,1)=str2double(tokens{2});
        SVstruct.sv_id(gc,1)=tokens(3);
        SVstruct.uid(gc,1)=strcat(sampleID,tokens(3));
        tokens= regexp(line,'(\w*)(\[|\])(\w*):(\d*)(\[|\])(\w*)','tokens','once');
	  if isempty(tokens),
	      continue
	  end
        SVstruct.altchr(gc,1)=chr_str2num(tokens(3));
        SVstruct.altpos(gc,1)=str2double(tokens{4});
        if isempty(tokens{1})
            SVstruct.strand(gc,1)='-';
        else
            SVstruct.strand(gc,1)='+';
        end
        if strcmp(tokens(5),']')
            SVstruct.altstrand(gc,1)='+';
        else
            SVstruct.altstrand(gc,1)='-';
        end
        SVstruct.dcc_project_code{gc,1}='DIPG';
        SVstruct.homseq{gc,1}=char(regexp(line,'(?<=\;HOMSEQ\=)\w+','match'));
        SVstruct.homlen(gc,1)=length(SVstruct.homseq{gc,1});           
        SVstruct.insertion{gc,1}=char(regexp(line,'(?<=\;INSERTION\=)\w+','match'));
        SVstruct.inslen(gc,1)=length(SVstruct.insertion{gc,1});
        %SVstruct.contig{gc,1}=char(regexp(line,'(?<=\;SCTG\=)\w+','match'));
        SVstruct.contig{gc,1}=char(regexp(line,'(?<=\;SCTG\=)[^;]+','match'));
        gc=gc+1;
    end
end

%[usid,ia_sid,ic_sid] = unique(SVstruct.sv_id,'stable');
%SVstruct_dd=SVstruct(ia_sid,:);
SVstruct_dd = SVstruct;


