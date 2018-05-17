class IdNotInSquareException(Exception):

    def __init__(self, square, sqid):
        self.square = square
        self.sqid = sqid

    def __str__(self):
        return "There is no "+str(self.square)+" cell in Square "+str(self.sqid)

