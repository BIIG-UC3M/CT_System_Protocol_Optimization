function d=readSimpleBinSliGap(fName,sizeRow,sizeCol,iniSlice,sizeSlices,fmt,headerSize,gapSize)

d = zeros([sizeRow,sizeCol,sizeSlices]);
sizeSli = sizeRow*sizeCol;

f=fopen(fName,'rb'); 

if(exist('headerSize','var'))
    fseek(f,headerSize,'bof');
end;

if(exist('fmt','var'))
    switch lower(fmt)
        case {'int8','uint8'}
            offset = (iniSlice-1)*sizeRow*sizeCol;
        case {'int16','uint16'}
            offset = (iniSlice-1)*sizeRow*sizeCol*2;
        case {'int32','uint32','single','float32'}
            offset = (iniSlice-1)*sizeRow*sizeCol*4;
        case {'int64','uint64','double','float64'}
            offset = (iniSlice-1)*sizeRow*sizeCol*8;
        otherwise
            offset = 0;
    end
    fseek(f,offset,'cof');
    for nSli = 1:sizeSlices,
        ptr = ((nSli-1)*sizeSli)+1;
        d(ptr:ptr+sizeSli-1)=fread(f,sizeSli,fmt);
        fseek(f,gapSize,'cof');
    end
else
   offset = (iniSlice-1)*sizeRow*sizeCol*4;
   fseek(f,offset,'cof');
    for nSli = 1:sizeSlices,
        fseek(f,gapSize,'cof');
        ptr = ((nSli-1)*sizeSli)+1;
        d(ptr:ptr+sizeSli-1)=fread(f,sizeSli,fmt);
    end
end;

d=reshape(d,[sizeRow sizeCol sizeSlices]);

fclose(f);