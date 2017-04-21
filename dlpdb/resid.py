
class ResID:
    def __init__(self, setChainID, setSeqNum, setICode):
        self.chainID = setChainID
        self.seqNum = setSeqNum
        self.iCode = setICode
    def __lt__(self, other):
        return ((self.chainID, self.seqNum, self.iCode) <
                (other.chainID, other.seqNum, other.iCode))
    def __le__(self, other):
        return not other < self
    def __gt__(self, other):
        return other < self
    def __ge__(self, other):
        return not self < other
    def __eq__(self, other):
        return (not self < other) and (not other < self)
    def __ne__(self, other):
        return (self < other) or (other < self)

    def __hash__(self):
        numChainIDs = 128
        numICodes = 128
        i = self.seqNum
        i *= numChainIDs
        i += ord(self.chainID)
        i *= numICodes
        i += ord(self.iCode)
        return i

    def __str__(self):
        return self.chainID+str(self.seqNum)+self.iCode

    def __repr__(self):
        return str(self)
