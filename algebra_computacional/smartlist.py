class SmartList(list):
    """docstring for SmartList"""
    def __init__(self, shape):
        super(SmartList, self).__init__()
        self.shape = shape
        self.list = None
        for i in xrange(len(shape)-1, -1, -1):
            self.list = [deepcopy(self.list) for x in range(shape[i])]
    def __getitem__(self, position):
        return self._getitemaux(self.list, position)
    def _getitemaux(self, l, position):
        if isinstance(position,int):
            return l[position]
        if len(position) > 1:
            return self._getitemaux(l[position[0]], position[1:])
        return l[position[0]]
    def __setitem__(self, position, value):
        return self._setitemaux(self.list, position, value)
    def _setitemaux(self, l, position, value):
        if isinstance(position,int):
            l[position] = value
            return
        if len(position) > 1:
            return self._setitemaux(l[position[0]], position[1:], value)
        l[position[0]] = value
    def __detitem__(self, position):
        return self._detitemaux(self.list, position)
    def _detitemaux(self, l, position):
        if isinstance(position,int):
            del l[position]
        if len(position) > 1:
            return self._detitemaux(l[position[0]], position[1:])
        print l
        del l[position[0]]