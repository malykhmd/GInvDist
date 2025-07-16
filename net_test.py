from net.ptp import P2P

server = P2P(60000, 2)
server.create_session('127.0.0.1')

server.send('127.0.0.1', 'test message')
