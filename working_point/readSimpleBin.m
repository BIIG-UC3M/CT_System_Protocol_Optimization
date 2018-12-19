function d=readSimpleBin(fName,sizeRow,sizeCol,sizeSlices,fmt,gapSize,headerSize)

f=fopen(fName,'rb');
d = zeros([sizeRow,sizeCol,sizeSlices]);
sizeSli = sizeRow*sizeCol;
if(exist('headerSize','var'))
    fseek(f,headerSize,'bof');
end;

if(exist('fmt','var'))
     for nSli = 1:sizeSlices,
         %keyboard
         %gapSize = cast(gapSize,'double')
        fseek(f,gapSize,'cof');
        ptr = ((nSli-1)*sizeSli)+1;
        d(ptr:ptr+sizeSli-1)=fread(f,sizeSli,fmt);
    end
    %d=fread(f,sizeRow*sizeCol*sizeSlices,fmt);
else
    d=fread(f,sizeRow*sizeCol*sizeSlices,'float32');

end;

d=reshape(d,[sizeRow sizeCol sizeSlices]);

fclose(f);