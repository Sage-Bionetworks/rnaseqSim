
import re

attribute_field = re.compile(r'(\w+) \"(.*)\"')


class Entry:
    def __init__(self, reference, source, feature, start, end, score, strand, frame, attribute) :
        self.reference = reference
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attribute = self.attribute_parse(attribute)

    def attribute_parse(self, attribute_line):
        tmp = attribute_line.split("; ")
        props = {}
        for i in tmp:
            res = attribute_field.search(i)
            if res is None:
                print "Can't parse: '%s'" % (i)
            else:
                props[res.group(1)] = res.group(2)
        return props

    def to_dict(self):
        return {
            "reference" : self.reference,
            "source" : self.source,
            "feature" : self.feature,
            "start" : self.start,
            "end" : self.end,
            "score" : self.score,
            "strand" : self.strand,
            "frame" : self.frame,
            "attribute" : self.attribute,
        }
    
    def __repr__(self):
        return str(self.to_dict())
    
            

class GTF:
    def __init__(self):
        self.entries = []
    
    def read(self, handle):
        for line in handle:
            if not line.startswith("#!"):
                tmp = line.rstrip().split("\t")
                entry = Entry(*tmp)
                if entry.source == "protein_coding":
                    self.entries.append( entry )

    def transcript_collect(self):
        out = {}
        for r in self.entries:
            if 'transcript_id' in r.attribute:
                k = r.attribute['transcript_id']
                if k not in out:
                    out[k] = []
                out[k].append(r)
        return out
            
    def __iter__(self):
        return self.entries.__iter__()
